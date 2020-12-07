import re


def paste_doc(lookup_dict):
    """
    Insert variable into obj's docstring.
    use `${val}` to insert a variable
    """

    def inner(obj):
        old_doc = obj.__doc__
        doc = ""
        old_e = 0
        for match in re.finditer("\${.*?}", old_doc):
            s, e = match.start(), match.end()
            doc += old_doc[old_e:s]
            content = old_doc[s + 2:e - 1]
            val = lookup_dict.get(content, "")
            doc += val
            old_e = e
        doc += old_doc[old_e:]
        obj.__doc__ = doc
        return obj

    return inner
