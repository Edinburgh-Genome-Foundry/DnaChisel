import hashlib
import base64

def install_extras_message(libname):
    return (
        "Could not load %s (is it installed ?). You can install it separately "
        " with:  pip install %s\n\n"
        "Install all dependencies for generating DNA Chisel reports with:"
        "\n\npip install dnachisel[reports]"
        % (libname, libname.lower().replace(" ", "_"))
    )

def file_hash(path=None, content=None, chars=7):
    if path is not None:
        with open(path, 'r') as f:
            content = f.read()
    hasher = hashlib.md5(content)
    md5bytes = hasher.digest()
    result = base64.urlsafe_b64encode(md5bytes).decode('ascii')
    return result[:chars]