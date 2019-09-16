Genbank API
===========

This section documents the Genbank API for all builtin specifications of
DNA chisel. 

.. contents::

Pattern removal
---------------

To remove a pattern in a given region, use the ``@no()`` label
(alias for ``@AvoidPattern()``). For instance to remove any GC pattern:

.. raw:: html

    <img class='annotation-example'
    src='../_static/images/genbank_annotations/avoid_pattern_cg.png'></img>

This method also supports IUPAC nucleotide mutation for degenerate sequences
(NKY etc.) For instance to remove any ATC and ATT (which could be written
together as ATH):

.. raw:: html

    <img class='annotation-example'
    src='../_static/images/genbank_annotations/avoid_pattern_noisoleucine.png'></img>

It is also possible to provide an enzyme restriction site by suffixing the enzyme
name with ``_site``:

.. raw:: html

    <img class='annotation-example'
    src='../_static/images/genbank_annotations/avoid_pattern_enzyme.png'></img>

Other pattern shorthands can be used for instance to find sequence repeats (here
we look for dimers repeated four times in a row):

.. raw:: html

    <img class='annotation-example'
    src='../_static/images/genbank_annotations/avoid_pattern_kmer.png'></img>

When no pattern shorthand cuts it, use a regular expression! The annotation
below forbids any TAA sequence following an ATG sequence.


.. raw:: html

    <img class='annotation-example'
    src='../_static/images/genbank_annotations/avoid_pattern_regex.png'></img>

Finally, note that you can also use ``no()`` as an objective by prefixing it
with a tilde ``~no()``, at which case the pattern may not be completely
eliminated, but its number of occurences will be minimized by the algorithm:

.. raw:: html

    <img class='annotation-example'
    src='../_static/images/genbank_annotations/avoid_pattern_cg_obj.png'></img>

Sequence protection
-------------------

Pattern insertion
-----------------



CDS and Codon optimization
--------------------------

GC content
----------

