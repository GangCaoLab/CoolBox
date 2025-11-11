import abc
import typing as T
import subprocess as subp
from os import path as osp

import numpy as np
import oxbow as ox
import pandas as pd

from ..logtools import get_logger
from ..filetool import opener, to_string
from ..genome import GenomeRange
from ..cmd import ensure_unix, ensure_tool_installed

log = get_logger(__name__)

BED12_FIELDS = ["chrom", "start", "end", "name", "score", "strand", "thick_start", "thick_end", "rgb", "block_count", "block_sizes", "block_starts"]

FMT2COLUMNS = {
    "bed6": BED12_FIELDS[:6],
    "bed9": BED12_FIELDS[:9],
    "bed12": BED12_FIELDS,
    "gtf": ["seqname", "source", "type", "start", "end", "score", "strand", "frame", "attributes"],
    "bigwig": ["chrom", "start", "end", "value"],
    "bedgraph": ["chrom", "start", "end", "value"],
    "bam": ["qname", "flag", "rname", "pos", "mapq", "cigar", "rnext", "pnext", "tlen", "seq", "qual", "tags"],
    "pairs": ["name", "chr1", "pos1", "chr2", "pos2", "strand1", "strand2"],
    "bedpe": ["chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score", "strand1", "strand2"],
}


def tabix_query(bgz_file, query: GenomeRange, split=True):
    """Call tabix and generate an array of strings for each line it returns."""
    ensure_unix()
    ensure_tool_installed("tabix")
    p = subp.Popen(['tabix', '-f', bgz_file, str(query)], stdout=subp.PIPE)
    for line in p.stdout:
        line = line.decode('utf-8')
        if not split:
            yield line
        else:
            yield line.strip().split('\t')


def pairix_query(bgz_file, query: GenomeRange, second: T.Optional[GenomeRange] = None,
                 open_region: bool = False, split: bool = True):
    ensure_unix()
    ensure_tool_installed("pairix")
    if second:
        query = f"{query}|{second}"
    else:
        if open_region:
            query = f"{query}|{query.chrom}"
    cmd = ['pairix', str(bgz_file), str(query)]
    p = subp.Popen(cmd, stdout=subp.PIPE)
    for line in p.stdout:
        line = line.decode('utf-8')
        if not split:
            yield line
        else:
            yield line.strip().split('\t')


def _is_bam_sorted(bam_path):
    p = subp.Popen(['samtools', 'view', '-H', bam_path], stdout=subp.PIPE)
    for line in p.stdout:
        line = line.decode("utf-8")
        if "SO:unsorted" in line:
            return False
    return True


def process_bam(bam_path):
    """Sort and index a BAM file.
    If input is a SAM file, it will be converted to a BAM file first."""
    if bam_path.endswith(".bam"):
        bai_path = bam_path + '.bai'
        if osp.exists(bai_path):
            return bam_path
        if osp.exists(bam_path[:-4] + '.sorted.bam.bai'):
            return bam_path[:-4] + '.sorted.bam.bai'

        ensure_unix()
        ensure_tool_installed("samtools")
        if not _is_bam_sorted(bam_path):
            sorted_bam_path = bam_path[:-4] + '.sorted.bam'
            subp.check_call(['samtools', 'sort', bam_path, '-o', sorted_bam_path])
        else:
            sorted_bam_path = bam_path
        subp.check_call(['samtools', 'index', sorted_bam_path])
    elif bam_path.endswith(".sam"):
        ensure_unix()
        ensure_tool_installed("samtools")
        sorted_bam_path = bam_path[:-4] + '.sorted.bam'
        subp.check_call(['samtools', 'sort', bam_path, '-o', sorted_bam_path])
        subp.check_call(['samtools', 'index', sorted_bam_path])
    else:
        raise IOError("BAM input file should be in .bam or .sam format")
    return sorted_bam_path


def query_bam(bam_path: str, gr: GenomeRange, split: bool = True):
    """Call tabix and generate an array of strings for each line it returns."""
    ensure_unix()
    ensure_tool_installed("samtools")
    p = subp.Popen(['samtools', 'view', bam_path, str(gr)], stdout=subp.PIPE)
    for line in p.stdout:
        line = line.decode('utf-8')
        if not split:
            yield line
        else:
            items = line.strip().split('\t')
            items_ = items[:11] + ["\t".join(items[12:])]
            items_[1] = int(items_[1])  # flag
            items_[3] = int(items_[3])  # pos
            items_[4] = int(items_[4])  # mapq
            yield items_


def _parse_samtools_cov(lines):
    covs = {}
    for line in lines[1:-1]:
        left, mid, _ = line.split("│")
        percent = float(left.strip("> %"))
        for i, c in enumerate(mid):
            covs.setdefault(i, 0)
            if c != ' ' and covs[i] == 0:
                covs[i] = percent
    covs = [covs[i] for i in sorted(covs.keys())]
    return covs


def coverage_by_samtools(bam_path, region, bins):
    ensure_unix()
    ensure_tool_installed("samtools")
    cmd = ["samtools", "coverage", bam_path, "-r", region, "-w", str(bins)]
    p = subp.Popen(cmd, stdout=subp.PIPE)
    lines = []
    for line in p.stdout:
        line = line.decode('utf-8')
        lines.append(line)
    covs = _parse_samtools_cov(lines)
    return np.array(covs)


def guess_bed_type(path: str) -> str:
    def get_no_comment_line(iter):
        line = next(iter)
        line = to_string(line)
        if line.startswith("#") or line.startswith("track") or \
                line.startswith("browser") or line.strip() == '':
            line = get_no_comment_line(iter)
        return line

    with opener(path) as f:
        try:
            fields = get_no_comment_line(iter=f)
        except StopIteration:
            raise ValueError(f"File is empty: {path}")
        fields = to_string(fields)
        line_values = fields.split("\t")

        if len(line_values) == 3:
            file_type = 'bed3'
        elif len(line_values) == 4:
            file_type = 'bedgraph'
        elif len(line_values) == 6:
            file_type = 'bed6'
        elif len(line_values) == 12:
            file_type = 'bed12'
        elif len(line_values) == 9:
            # this is a case where a specific color is encoded in the 10 field of the bed file
            file_type = 'bed9'
        elif len(line_values) > 6:
            # assume bed6
            file_type = 'bed6'
            log.warning("Number of fields in BED file is not standard. Assuming bed6\n")
        else:
            raise ValueError(f"Number of fields in BED file is not standard: {len(line_values)}")
    return file_type


def get_columns(path: str) -> T.List[str]:
    """Get the columns for the file."""
    if path.endswith(".bgz"):
        _p = path.rstrip(".bgz")
    else:
        _p = path
    suffix = osp.splitext(_p)[1].lower()
    if suffix in [".bed", ".bedgraph", ".bg"]:
        return FMT2COLUMNS[guess_bed_type(path)]
    else:
        fmt = suffix[1:]
        if fmt == 'bw':
            fmt = 'bigwig'
        if fmt == 'bg':
            fmt = 'bedgraph'
        if fmt in FMT2COLUMNS:
            return FMT2COLUMNS[fmt]
        else:
            raise ValueError(f"There is no predefined columns for file type: {fmt}")


class TabFileReader(abc.ABC):
    """Generic tab-separated file reader.

    Including:
    - BED
    - BedGraph
    - bigWig
    - bigBED
    - GTF
    - BAM
    """
    def __init__(self, path: str, columns: T.List[str] | None = None, **params):
        self.path = path
        if path.endswith(".bgz"):
            _p = path.rstrip(".bgz")
        else:
            _p = path
        suffix = osp.splitext(_p)[1].lower()
        self.suffix = suffix
        if columns is None:
            columns = get_columns(path)
        self.columns = columns
        self.params = params
        self.is_2d = False
        if suffix in [".bedpe", ".pairs"]:
            self.is_2d = True

    @abc.abstractmethod
    def query(self, gr: GenomeRange, **kwargs) -> pd.DataFrame:
        """Query the file"""
        pass

    def query_var_chr(self, gr: GenomeRange, **kwargs) -> pd.DataFrame:
        """Query the file with variable chromosome names.
        
        First try to query with the original chromosome names, if no results,
        try to query with the variable chromosome names."""
        df = self.query(gr, **kwargs)
        if df.shape[0] > 0:
            return df
        else:
            gr = gr.change_chrom_names()
            gr2 = kwargs.get("second", None)
            if gr2:
                gr2 = gr2.change_chrom_names()
                kwargs["second"] = gr2
            df = self.query(gr, **kwargs)
            return df


def _dict_to_gtf_attr(attr_dict):
    """Convert a dictionary to a GTF attribute string. 
    For recovery the GTF attribute column fetched by TabFileReaderWithOxbow. """
    valid_items = {
        k: v for k, v in attr_dict.items()
        if v is not None and v != ''
    }
    return ' '.join([f'{k} "{v}";' for k, v in valid_items.items()])


class TabFileReaderWithOxbow(TabFileReader):
    def __init__(self, path: str, columns: T.List[str] | None = None, **params):
        super().__init__(path, columns, **params)
        suffix = self.suffix
        if suffix == ".gtf":
            ds = ox.from_gtf(self.path)
        elif suffix in [".bed", ".bedgraph", ".bg"]:
            ds = ox.from_bed(self.path)
        elif suffix in ['.bw', '.bigwig']:
            ds = ox.from_bigwig(self.path)
        elif suffix == '.bam':
            ds = ox.from_bam(self.path)
        else:
            raise NotImplementedError(f"Unsupported file type: {suffix}")
        self.ds = ds

    def query(self, gr: GenomeRange, **kwargs) -> pd.DataFrame:
        sub = self.ds.regions(str(gr))
        try:
            df = sub.pd()
            if self.suffix in [".bed", ".bedgraph", ".bg"]:
                rest = df.pop('rest')
                df_rest = rest.str.split('\t', expand=True)
                df = pd.concat([df, df_rest], axis=1)
                df.columns = self.columns
            elif self.suffix == ".bam":
                if 'end' in df.columns:
                    df.pop('end')
            elif self.suffix == ".gtf":
                if 'seqid' in df.columns:
                    df.rename(columns={'seqid': 'seqname'}, inplace=True)
                if 'attributes' in df.columns:
                    df['attributes'] = df['attributes'].apply(_dict_to_gtf_attr)
            df = _convert_dtype(df)
        except ValueError as e:
            # empty region
            log.error(str(e))
            df = pd.DataFrame(columns=self.columns)
        return df


class TabFileReaderWithTabix(TabFileReader):
    def __init__(self, path: str, columns: T.List[str] | None = None, **params):
        super().__init__(path, columns, **params)
        suffix = self.suffix
        ensure_unix()
        if suffix in [".bedpe", ".pairs"]:
            ensure_tool_installed("pairix")
        else:
            ensure_tool_installed("tabix")

    def query(
            self,
            gr: GenomeRange,
            second: T.Optional[GenomeRange] = None,
            open_region: bool = False,
            **kwargs) -> pd.DataFrame:
        if self.is_2d:
            itr = pairix_query(self.path, gr, second, open_region, split=True)
        elif self.suffix == '.bam':
            itr = query_bam(self.path, gr, split=True)
        else:
            itr = tabix_query(self.path, gr, split=True)
        rows = list(itr)
        if len(rows) > 0:
            n_cols = len(rows[0])
            columns = self.columns[:n_cols]
        else:
            columns = self.columns
        df = pd.DataFrame(rows, columns=columns)
        df = _convert_dtype(df)
        return df


def _convert_dtype(df: pd.DataFrame) -> pd.DataFrame:
    def _convert_tp(col_name: str, dtype):
        try:
            df[col_name] = df[col_name].astype(dtype)
        except ValueError:
            pass

    # convert non-string chromosome columns to string
    for chr_col_name in ['chrom', 'seqname', 'rname', 'chr1', 'chr2', 'chrom1', 'chrom2']:
        if chr_col_name in df.columns:
            dtype = df[chr_col_name].dtype
            if dtype != 'object':
                _convert_tp(chr_col_name, str)
    # convert integer columns to int
    for col_name in ['start', 'end', 'pos1', 'pos2', 'start1', 'start2', 'end1', 'end2', 'block_count', 'thick_start', 'thick_end']:
        if col_name in df.columns:
            dtype = df[col_name].dtype
            if dtype != int:
                _convert_tp(col_name, int)
    # convert float columns to float
    for col_name in ['value', 'score']:
        if col_name in df.columns:
            dtype = df[col_name].dtype
            if dtype != float:
                _convert_tp(col_name, float)
    return df


class TabFileReaderInMemory(TabFileReader):
    def __init__(self, path: str, columns: T.List[str] | None = None, **params):
        super().__init__(path, columns, **params)
        with opener(path) as f:
            self.df = pd.read_csv(f, sep='\t', comment='#')
            self.df.columns = self.columns[:len(self.df.columns)]
        self.df = _convert_dtype(self.df)

    def query(
            self,
            gr: GenomeRange,
            second: T.Optional[GenomeRange] = None,
            open_region: bool = False,
            **kwargs) -> pd.DataFrame:
        if self.is_2d:
            if self.suffix == '.pairs':
                field_names = {
                    'chr1': 'chr1',
                    'start1': 'pos1',
                    'end1': 'pos1',
                    'chr2': 'chr2',
                    'start2': 'pos2',
                    'end2': 'pos2',
                }
            else:
                # bedpe
                field_names = {
                    'chr1': 'chrom1',
                    'start1': 'start1',
                    'end1': 'end1',
                    'chr2': 'chrom2',
                    'start2': 'start2',
                    'end2': 'end2',
                }

            if second is not None:
                sdf = self.df.query(
                    f"{field_names['chr1']} == '{gr.chrom}' and {field_names['start1']} >= {gr.start} and {field_names['end1']} <= {gr.end} "
                    f"and {field_names['chr2']} == '{second.chrom}' and {field_names['start2']} >= {second.start} and {field_names['end2']} <= {second.end}"
                )
                return sdf
            else:
                if open_region:
                    q = (
                        f"{field_names['chr1']} == '{gr.chrom}' and {field_names['start1']} >= {gr.start} and {field_names['end1']} <= {gr.end} "
                        f"and {field_names['chr2']} == '{gr.chrom}'"
                    )
                    sdf = self.df.query(q)
                else:
                    sdf = self.df.query(
                        f"{field_names['chr1']} == '{gr.chrom}' and {field_names['start1']} >= {gr.start} and {field_names['end1']} <= {gr.end} "
                        f"and {field_names['chr2']} == '{gr.chrom}' and {field_names['start2']} >= {gr.start} and {field_names['end2']} <= {gr.end}"
                    )
                return sdf
        else:
            chrom_str = 'chrom'
            if self.suffix == '.gtf':
                chrom_str = 'seqname'
            sdf = self.df.query(f"{chrom_str} == '{gr.chrom}' and start >= {gr.start} and end <= {gr.end}")
            return sdf


def _build_bgz_file(
        path: str,
        col_chrom: T.Optional[int] = None,
        col_start: T.Optional[int] = None,
        col_end: T.Optional[int] = None) -> str:
    input_is_gz = False
    if path.endswith(".bgz"):
        prefix = path[:-4]
    elif path.endswith(".gz"):
        input_is_gz = True
        prefix = path[:-3]
    else:
        prefix = path
    output_path = prefix + ".bgz"
    if osp.exists(output_path):
        log.info(f"Bgz file already exists, skip building:\n{output_path}")
        return output_path

    ensure_unix()
    ensure_tool_installed("bgzip")

    cat_cmd = "zcat" if input_is_gz else "cat"
    if prefix.lower().endswith(".gtf"):
        cmd = f'{cat_cmd} {path} | grep -v ^"#" | sort -k1,1 -k4,4n | bgzip > {output_path}'
    elif prefix.lower().endswith('.bed') or prefix.lower().endswith('.bedgraph') or prefix.lower().endswith('.bg'):
        cmd = f'{cat_cmd} {path} | sort -k1,1 -k2,2n | bgzip > {output_path}'
    elif prefix.lower().endswith('.bedpe'):
        cmd = f'{cat_cmd} {path} | sort -k1,1 -k4,4 -k2,2n -k5,5n | bgzip > {output_path}'
    elif prefix.lower().endswith('.pairs'):
        cmd = f'{cat_cmd} {path} | grep -v ^"#" | sort -k2,2 -k4,4 -k3,3n -k5,5n | bgzip > {output_path}'
    else:
        assert col_chrom is not None, "col_chrom is required"
        assert col_start is not None, "col_start is required"
        c_c, c_s = col_chrom + 1, col_start + 1
        if col_end is None:
            cmd = f'{cat_cmd} {path} | grep -v ^"#" | sort -k{c_c},{c_c} -k{c_s},{c_s}n | bgzip > {output_path}'
        else:
            c_e = col_end + 1
            cmd = f'{cat_cmd} {path} | grep -v ^"#" | sort -k{c_c},{c_c} -k{c_s},{c_s}n -k{c_e},{c_e}n | bgzip > {output_path}'
    log.info(f"Build bgz file, save to {output_path}")
    log.info(f"Command:\n{cmd}")
    subp.check_call(cmd, shell=True)
    return output_path


def _index_bgz_file(
        bgz_path: str,
        col_chrom: T.Optional[int] = None,
        col_start: T.Optional[int] = None,
        col_end: T.Optional[int] = None) -> str:
    index_file = bgz_path + ".tbi"
    if osp.exists(index_file):
        log.info(f"Tabix index already exists, skip indexing:\n{index_file}")
        return index_file

    ensure_unix()
    ensure_tool_installed("tabix")
    prefix = bgz_path[:-4]
    if prefix.lower().endswith(".gtf"):
        cmd = ['tabix', '-p', 'gff', bgz_path]
    elif prefix.lower().endswith('.bed'):
        cmd = ['tabix', '-0', '-p', 'bed', bgz_path]
    elif prefix.lower().endswith('.bedgraph') or prefix.lower().endswith('.bg'):
        cmd = ['tabix', '-0', '-b', '2', '-e', '3', bgz_path]
    elif prefix.lower().endswith('.bedpe'):
        ensure_tool_installed("pairix")
        cmd = ['pairix', '-f', '-s', '1', '-d', '4', '-b', '2', '-e', '3', '-u', '5', '-v', '6', bgz_path]
    elif prefix.lower().endswith('.pairs'):
        ensure_tool_installed("pairix")
        cmd = ['pairix', '-f', '-p', 'pairs', bgz_path]
    else:
        assert col_chrom is not None, "col_chrom is required"
        assert col_start is not None, "col_start is required"
        c_c, c_s = col_chrom + 1, col_start + 1
        if col_end is None:
            cmd = ['tabix', '-0', '-s', str(c_c), '-b', str(c_s), '-e', str(c_s), bgz_path]
        else:
            c_e = col_end + 1
            cmd = ['tabix', '-0', '-s', str(c_c), '-b', str(c_s), '-e', str(c_e), bgz_path]
    log.info(f"Index bgz file, save to {index_file}")
    log.info(f"Command:\n{' '.join(cmd)}")
    subp.check_call(cmd)
    return index_file


def index_tab_file(
        path: str,
        col_chrom: T.Optional[int] = None,
        col_start: T.Optional[int] = None,
        col_end: T.Optional[int] = None) -> str:
    if path.endswith(".bam"):
        return process_bam(path)
    elif path.endswith(".bw") or path.lower().endswith(".bigwig"):
        return path
    else:
        bgz_path = _build_bgz_file(path, col_chrom, col_start, col_end)
        _index_bgz_file(bgz_path, col_chrom, col_start, col_end)
        return bgz_path


def get_indexed_tab_reader(
        path: str,
        columns: T.Optional[T.List[str]] = None) -> TabFileReader:
    if columns is None:
        try:
            columns = get_columns(path)
        except ValueError:
            raise ValueError(f"You must specify the columns for this file: {path}")
    col_chrom = columns.index("chrom") if "chrom" in columns else None
    if 'start' in columns:
        col_start = columns.index("start")
    elif 'pos' in columns:
        col_start = columns.index("pos")
    else:
        col_start = None
    col_end = columns.index("end") if "end" in columns else None
    try:
        indexed_path = index_tab_file(path, col_chrom, col_start, col_end)
        reader = TabFileReaderWithOxbow(indexed_path, columns)
        return reader
    except OSError as e:
        log.error(str(e))
        log.warning(f"Try to use TabFileReaderInMemory instead")
        reader = TabFileReaderInMemory(path, columns=columns)
        return reader
    except NotImplementedError:
        # Unsupported file type for oxbow
        try:
            reader = TabFileReaderWithTabix(indexed_path, columns=columns)
            return reader
        except OSError as e:
            log.error(str(e))
            log.warning(f"Try to use TabFileReaderInMemory instead")
            reader = TabFileReaderInMemory(path, columns=columns)
            return reader

