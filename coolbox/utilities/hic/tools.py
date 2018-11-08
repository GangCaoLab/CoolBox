
def infer_resolution(genome_range, resolutions, bin_thresh=1000):
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

