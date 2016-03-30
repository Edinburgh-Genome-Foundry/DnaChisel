


DNAChisel
===========


.. raw:: html

    <a href="https://twitter.com/share" class="twitter-share-button"
    data-text="DNAChisel - A Python module for printing with living matter" data-size="large" data-hashtags="Bioprinting">Tweet
    </a>
    <script>!function(d,s,id){var js,fjs=d.getElementsByTagName(s)[0],p=/^http:/.test(d.location)?'http':'https';
    if(!d.getElementById(id)){js=d.createElement(s);js.id=id;js.src=p+'://platform.twitter.com/widgets.js';
    fjs.parentNode.insertBefore(js,fjs);}}(document, 'script', 'twitter-wjs');
    </script>
    <iframe src="http://ghbtns.com/github-btn.html?user=Edinburgh-Genome-Foundry&repo=dnachisel&type=watch&count=true&size=large"
    allowtransparency="true" frameborder="0" scrolling="0" width="152px" height="30px" margin-bottom="30px"></iframe>

DNAChisel is a Python module to produce living art. It transforms an image into files that a liquid dispenser can use to *print* the image to a plate using pigmented yeast or bacteria.

Here are two examples of bio-art:

.. figure:: images/bioprint_dolly.jpeg
    :align: center

    Dolly drawn with baker yeast (white), violacein-producing yeast (black), and carotene-producing yeast (orange)


.. figure:: images/bioprint_england.jpeg
    :align: center

    England flag drawn with 3 different strains of the bacterium *E. coli*.


DNAChisel is released on Github_ under the MIT licence (¢ Edinburgh Genome Foundry), everyone is welcome to contribute !

DNAChisel was written at the Edinburgh Genome Foundry by Zulko_ after an original idea and Matlab code by Mike Shen (`Mike's project on Github <https://github.com/mshen5/BioPointillism>`_).




Installation
--------------

If you have PIP installed: ::

    (sudo) pip install ez_setup dnachisel

Or unzip the source code in a directory and type in a terminal: ::

    sudo python setup.py install


Usage
--------

In the same folder as your code, place an image. It can have any resolution, but keep in mind that the width/height ratio of the plate it is printed on is 1.5. Make sure that a specific color is used to mark the un-pigmented background of the image, here we use blue:

.. figure:: images/dolly.jpeg
    :align: center

Then write the following code in ``dolly.py``:
::
    from DNAChisel import bioprint

    bioprint(
        image_filename="../docs/images/dolly.jpeg",
        output_filename="dolly.csv",
        bg_color=[0, 0, 255], # blue background represents empty wells
        pigments_wells={"A1": [0, 0, 0],  # black yeast in source well A1
                        "A2": [250, 120, 10],  # orange yeast in well A2
                        "A3": [255, 255, 255]},  # white yeast in well A3
        quantified_image_filename="dolly_preview.jpeg"
    )

Execute in a terminal with ``python dolly.py``. This will produce a ``dolly.csv`` file as well as a preview image of the final printing (so that you can check if the image looks good at this low resolution).

.. figure:: images/dolly_preview.png
    :align: center

Prepare a source plate with the right pigmented yeasts in wells A1, A2, A3, use an agar plate as the destination plate, and feed ``dolly.csv`` to the `Labcyte Echo <http://www.labcyte.com/products/liquidhandling/echo-555-liquid-handler>`_. Once the printing is finished, incubate 2 days at 30C (it would be one day at 37C for bacteria). Enjoy the result !


Reference manual
-----------------

.. autofunction:: DNAChisel.bioprint

.. _Zulko: https://github.com/Zulko/
.. _Github: https://github.com/Edinburgh-Genome-Foundry/DNAChisel
