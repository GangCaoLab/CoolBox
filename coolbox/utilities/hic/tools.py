from ..genome import GenomeRange


def hicmat_filetype(path):
    if path.endswith(".hic"):
        return '.hic'
    else:
        p = path.split("::")[0]
        if p.endswith((".cool", ".mcool")):
            return '.cool'
        else:
            raise NotImplementedError("File type of {} not supported for HicMat".format(path))


def infer_resolution(genome_range: GenomeRange, resolutions: object, bin_thresh: object = 1000) -> object:
    """
    Inference appropriate resolution.
    """
    resolutions.sort()
    reso = resolutions[0]
    for r in resolutions:
        num_bins = genome_range.length // reso
        if num_bins >= bin_thresh:
            reso = r
        else:
            break
    return reso


def is_multi_cool(cooler_file):
    """
    Judge a cooler is muliti-resolution cool or not.

    Parameters
    ----------
    cooler_file : str
        Path to cooler file.
    """
    import re
    if re.match(".+::.+$", cooler_file):
        return False

    import h5py
    h5_file = h5py.File(cooler_file, 'r')
    is_multi = 'pixels' not in h5_file  # use "pixels" group distinguish is multi-cool or not
    h5_file.close()
    return is_multi


def get_cooler_resolutions(cooler_file, is_multi=True):
    """
    Get the resolutions of a muliti-cooler file

    Parameters
    ----------
    cooler_file : str
        Path to cooler file.
    """
    import h5py
    h5_file = h5py.File(cooler_file, 'r')
    if is_multi:
        if 'resolutions' in h5_file:
            resolutions = list(h5_file['resolutions'])
            resolutions = [int(res) for res in resolutions]
        else:
            resolutions = [int(h5_file[i].attrs['bin-size']) for i in list(h5_file)]
        resolutions.sort()
        h5_file.close()
    else:
        resolutions = [int(h5_file.attrs['bin-size'])]
    return resolutions
