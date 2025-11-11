import platform
import pandas as pd
import pytest

from coolbox.utilities.reader.tab import (
    TabFileReaderInMemory,
    TabFileReaderWithTabix,
    TabFileReaderWithOxbow,
    get_indexed_tab_reader,
    guess_bed_type,
    FMT2COLUMNS,
    index_tab_file,
)
from coolbox.utilities.cmd import check_tool


def test_guess_bed_type(data_dir, test_itv):
    assert guess_bed_type(f"{data_dir}/bed_{test_itv}.bed") == "bed12"
    assert guess_bed_type(f"{data_dir}/bed6_{test_itv}.bed") == "bed6"
    assert guess_bed_type(f"{data_dir}/bed9_{test_itv}.bed") == "bed9"
    assert guess_bed_type(f"{data_dir}/bedgraph_{test_itv}.bg") == "bedgraph"


def test_inmemory_bed6_query(data_dir, test_interval, empty_interval, test_itv):
    path = f"{data_dir}/bed6_{test_itv}.bed"
    rdr = TabFileReaderInMemory(path)
    df = rdr.query(test_interval)
    assert isinstance(df, pd.DataFrame)
    assert list(df.columns) == FMT2COLUMNS["bed6"]
    assert (df["chrom"] == test_interval.chrom).all()
    assert (df["start"] >= test_interval.start).all()
    assert (df["end"] <= test_interval.end).all()
    # empty interval returns empty frame
    df_empty = rdr.query(empty_interval)
    assert df_empty.shape[0] == 0


def test_inmemory_bed12_query(data_dir, test_interval, test_itv):
    path = f"{data_dir}/bed_{test_itv}.bed"
    rdr = TabFileReaderInMemory(path)
    df = rdr.query(test_interval)
    assert isinstance(df, pd.DataFrame)
    assert list(df.columns) == FMT2COLUMNS["bed12"]
    assert df.shape[0] > 0


def test_inmemory_bedgraph_query(data_dir, test_interval, empty_interval, test_itv):
    path = f"{data_dir}/bedgraph_{test_itv}.bg"
    rdr = TabFileReaderInMemory(path)
    df = rdr.query(test_interval)
    assert list(df.columns) == FMT2COLUMNS["bedgraph"]
    assert df.shape[0] > 0
    rdr.query(empty_interval)  # should not raise


def test_inmemory_gtf_query(data_dir, test_interval, empty_interval, test_itv):
    path = f"{data_dir}/gtf_{test_itv}.gtf"
    rdr = TabFileReaderInMemory(path)
    df = rdr.query(test_interval)
    # GTF uses seqname instead of chrom in in-memory reader
    assert list(df.columns) == FMT2COLUMNS["gtf"]
    assert df.shape[0] > 0
    df_empty = rdr.query(empty_interval)
    assert df_empty.shape[0] == 0


def test_inmemory_pairs_query(data_dir, test_interval, test_itv):
    path = f"{data_dir}/pairs_{test_itv}.pairs"
    rdr = TabFileReaderInMemory(path)
    # 1D within same region
    df_same = rdr.query(test_interval)
    assert list(df_same.columns) == FMT2COLUMNS["pairs"]
    assert df_same.shape[0] > 0
    # 2D with explicit second region
    df_2d = rdr.query(test_interval, second=test_interval)
    assert df_2d.shape[0] > 0
    # open_region allows second chrom-only matching
    df_open = rdr.query(test_interval, open_region=True)
    assert df_open.shape[0] >= df_same.shape[0]
    assert df_open['pos1'].dtype == int
    assert df_open['pos2'].dtype == int


def test_inmemory_bedpe_query(data_dir, test_interval, test_itv):
    path = f"{data_dir}/bedpe_{test_itv}.bedpe"
    rdr = TabFileReaderInMemory(path)
    df_same = rdr.query_var_chr(test_interval)
    assert list(df_same.columns) == FMT2COLUMNS["bedpe"]
    assert df_same.shape[0] > 0
    df_2d = rdr.query_var_chr(test_interval, second=test_interval)
    assert df_2d.shape[0] > 0
    assert df_2d['start1'].dtype == int
    assert df_2d['start2'].dtype == int
    assert df_2d['end1'].dtype == int
    assert df_2d['end2'].dtype == int


@pytest.mark.skipif(platform.system() != "Windows", reason="Specific to Windows fallback behavior")
def test_get_indexed_tab_reader_fallback_to_inmemory_on_windows(data_dir, test_itv):
    # On Windows, bgzip/tabix are unavailable and ensure_unix raises, so it must fall back
    path = f"{data_dir}/bed6_{test_itv}.bed"
    rdr = get_indexed_tab_reader(path, columns=FMT2COLUMNS["bed6"])
    assert isinstance(rdr, TabFileReaderInMemory)

 
@pytest.mark.skipif(
    platform.system() == "Windows" or not check_tool("tabix")[0],
    reason="tabix not available or non-Unix",
)
def test_tabix_bedgraph_query(data_dir, test_interval):
    # Pre-indexed bedGraph bgz/tbi provided in test_data
    path = f"{data_dir}/chr9.1.pc.bedGraph"
    indexed_path = index_tab_file(path)
    rdr = TabFileReaderWithTabix(indexed_path)
    df = rdr.query(test_interval)
    assert list(df.columns) == FMT2COLUMNS["bedgraph"]
    assert df.shape[0] > 0
    assert df['start'].dtype == int
    assert df['end'].dtype == int
    assert df['value'].dtype == float


@pytest.mark.skipif(
    platform.system() == "Windows"
    or not (check_tool("bgzip")[0] and check_tool("pairix")[0]),
    reason="bgzip/pairix not available or non-Unix",
)
def test_tabix_pairs_query(data_dir, test_interval, test_itv):
    # Build .bgz and index via pairix, then query with Tabix reader
    path = f"{data_dir}/pairs_{test_itv}.pairs"
    bgz_path = index_tab_file(path)
    rdr = TabFileReaderWithTabix(bgz_path)
    df_same = rdr.query(test_interval)
    assert list(df_same.columns) == FMT2COLUMNS["pairs"]
    assert df_same.shape[0] > 0
    df_2d = rdr.query(test_interval, second=test_interval)
    assert df_2d.shape[0] > 0
    assert df_2d['pos1'].dtype == int
    assert df_2d['pos2'].dtype == int


def test_oxbow_bed_query(data_dir, test_interval, test_itv):
    # Oxbow-based reader for BED files
    path = f"{data_dir}/bed6_{test_itv}.bed"
    indexed_path = index_tab_file(path)
    rdr = TabFileReaderWithOxbow(indexed_path)
    df = rdr.query(test_interval)
    assert isinstance(df, pd.DataFrame)
    assert list(df.columns) == FMT2COLUMNS["bed6"]
    assert df.shape[0] > 0


def test_oxbow_gtf_query(data_dir, test_interval, test_itv):
    path = f"{data_dir}/gtf_{test_itv}.gtf"
    indexed_path = index_tab_file(path)
    rdr = TabFileReaderWithOxbow(indexed_path)
    df = rdr.query(test_interval)
    # ensure essential GTF columns exist after oxbow conversion
    for col in ["seqname", "start", "end", "strand"]:
        assert col in df.columns
    assert df.shape[0] > 0


def test_oxbow_bigwig_query(data_dir, test_interval, test_itv):
    path = f"{data_dir}/bigwig_{test_itv}.bw"
    rdr = TabFileReaderWithOxbow(path)
    df = rdr.query(test_interval)
    # BigWig should have chrom,start,end,value
    for col in FMT2COLUMNS["bigwig"]:
        assert col in df.columns
    assert df.shape[0] > 0


@pytest.mark.skipif(
    platform.system() == "Windows" or not check_tool("samtools")[0],
    reason="samtools not available or non-Unix",
)
def test_tabix_bam_query(data_dir, test_interval, test_itv):
    path = f"{data_dir}/bam_{test_itv}.bam"
    indexed_path = index_tab_file(path)
    rdr = TabFileReaderWithTabix(indexed_path)
    df = rdr.query(test_interval)
    assert list(df.columns) == FMT2COLUMNS["bam"]
    assert df.shape[0] > 0


@pytest.mark.skipif(
    platform.system() == "Windows" or not check_tool("tabix")[0],
    reason="tabix not available or non-Unix",
)
def test_tabix_gtf_query(data_dir, test_interval, test_itv):
    path = f"{data_dir}/gtf_{test_itv}.gtf"
    indexed_path = index_tab_file(path)
    rdr = TabFileReaderWithTabix(indexed_path)
    df = rdr.query(test_interval)
    assert list(df.columns) == FMT2COLUMNS["gtf"]
    assert df.shape[0] > 0


@pytest.mark.skipif(
    platform.system() == "Windows"
    or not (check_tool("bgzip")[0] and check_tool("pairix")[0]),
    reason="bgzip/pairix not available or non-Unix",
)
def test_tabix_bedpe_query(data_dir, test_interval, test_itv):
    path = f"{data_dir}/bedpe_{test_itv}.bedpe"
    indexed_path = index_tab_file(path)
    rdr = TabFileReaderWithTabix(indexed_path)
    df_same = rdr.query_var_chr(test_interval)
    assert list(df_same.columns) == FMT2COLUMNS["bedpe"]
    assert df_same.shape[0] > 0
    df_2d = rdr.query_var_chr(test_interval, second=test_interval)
    assert df_2d.shape[0] > 0
    assert df_2d['start1'].dtype == int
    assert df_2d['start2'].dtype == int
    assert df_2d['end1'].dtype == int
    assert df_2d['end2'].dtype == int


@pytest.mark.skipif(
    platform.system() == "Windows" or not check_tool("samtools")[0],
    reason="samtools not available or non-Unix",
)
def test_oxbow_bam_query(data_dir, test_interval, test_itv):
    path = f"{data_dir}/bam_{test_itv}.bam"
    indexed_path = index_tab_file(path)
    rdr = TabFileReaderWithOxbow(indexed_path)
    df = rdr.query(test_interval)
    assert df.shape[0] > 0
    # Check a few SAM/BAM fields commonly present
    for col in FMT2COLUMNS["bam"]:
        assert col in df.columns
