DNA Chisel
==========

DNAChisel is a Python library to modify the nucleotides of DNA sequences with respect to a set of
constraints and optimization objectives.


Installation
-------------

You can install DNAChisel through PIP
::
    sudo pip install dnachisel

Alternatively, you can unzip the sources in a folder and type
::
    sudo python setup.py install


Example
--------



.. code:: python

    from canvas = DNACanvas(
        sequence = "ATGCGTGTGCGTATGCGTGTGTGCGTGATG",
        constraints = [
            Pattern("ATTCTT", window = [100, 200]),
            NoPattern("AGTC", window = [300, 600]),

            PreserveORF(start1, end1),
            PreserveORF(start2, end2),
            PreserveORF(start3, end3),

            GCPercent(min=0.4, max=0.6, window=100)
            GCPercent(min=0.7, max=0.6, window=100)
        ]
        objective = [
            CodonOptimization(start=, end=, organism=)

        ]
    )
    max_local_gc_content(percentage=20, window=60)
    optimization= []


Contribute
----------

DNAChisel is an open-source library originally written at the
Edinburgh Genome Foundry by Zulko_.
It is released on Github under the MIT licence, everyone is welcome to contribute.
