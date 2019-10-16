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