def cm2inch(*tupl):
    """
    convert length unit from cm to inch.

    >>> cm2inch(10)
    3.937007874015748
    >>> cm2inch(2, 2)
    (0.7874015748031495, 0.7874015748031495)
    >>> cm2inch((1, 5))
    (0.39370078740157477, 1.968503937007874)
    """
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i / inch for i in tupl[0])
    else:
        if len(tupl) != 1:
            return tuple(i / inch for i in tupl)
        else:
            return tupl[0] / inch


def inch2cm(*tupl):
    """
    convert length unit from inch to cm.

    """
    cm = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i * cm for i in tupl[0])
    else:
        if len(tupl) != 1:
            return tuple(i * cm for i in tupl)
        else:
            return tupl[0] * cm


def rgb2hex(r, g, b):
    """
    Convert rgb color to hex format.

    >>> rgb2hex(129, 154, 70)
    '#819a46'
    >>> rgb2hex(-10, 256, -1)
    Traceback (most recent call last):
    ...
    AssertionError: (r, g, b) value must within range 0 ~ 255.
    """
    assert (0, 0, 0) <= (r, g, b) <= (255, 255, 255), \
        "(r, g, b) value must within range 0 ~ 255."
    hex = "#{:02x}{:02x}{:02x}".format(r, g, b)
    return hex


def hex2rgb(color_hex):
    """
    Convert hex color code to rgb tuple.

    >>> hex2rgb('#819a46')
    (129, 154, 70)
    """
    code = color_hex[1:]
    r = int(code[0:2], 16)
    g = int(code[2:4], 16)
    b = int(code[4:6], 16)
    return r, g, b


def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    This function is from: https://stackoverflow.com/a/20528097/8500469

    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax/(vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          `midpoint` and 1.0.

    '''
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt

    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False),
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap


def get_size(obj, seen=None):
    """
    Recursively finds size of objects

    From:
        https://stackoverflow.com/a/40880923/8500469
    """
    import sys
    size = sys.getsizeof(obj)
    if seen is None:
        seen = set()
    obj_id = id(obj)
    if obj_id in seen:
        return 0
    # Important mark as seen *before* entering recursion to gracefully handle
    # self-referential objects
    seen.add(obj_id)
    if isinstance(obj, dict):
        size += sum([get_size(v, seen) for v in obj.values()])
        size += sum([get_size(k, seen) for k in obj.keys()])
    elif hasattr(obj, '__dict__'):
        size += get_size(obj.__dict__, seen)
    elif hasattr(obj, '__iter__') and not isinstance(obj, (str, bytes, bytearray)):
        size += sum([get_size(i, seen) for i in obj])
    return size


def fig2bytes(fig, encode='svg', dpi=None):
    """
    Convert matplotlib.figure.Figure/IPython.display.SVG
    object to image bytes.
    """
    from IPython.display import SVG
    if isinstance(fig, SVG):
        img_bytes = fig.data.encode('utf-8')
    else:
        import io
        buf = io.BytesIO()
        fig.savefig(buf, format=encode, dpi=None)
        buf.seek(0)
        img_bytes = buf.read()
        buf.close()
    return img_bytes
