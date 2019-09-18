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

To prevent sections of the sequence to be modified, use ``@keep`` (alias for
``@AvoidChanges``):

.. raw:: html

    <img class='annotation-example'
    src='../_static/images/genbank_annotations/keep.png'></img>

You can also use ``keep`` as an optimization objective, at which case sequence
modifications will not be strictly forbidden, but they will be minimized:

.. raw:: html

    <img class='annotation-example'
    src='../_static/images/genbank_annotations/keep_obj.png'></img>

Pattern insertion
-----------------
You can control how many times a pattern should appear in a sequence region
with the ``@insert()`` specification (short form of ``@EnforcePatternOccurence``):

.. raw:: html

    <img class='annotation-example'
    src='../_static/images/genbank_annotations/insert.png'></img>

By default ``@insert()`` ensures that exactly one occurence of the pattern is
present in the given region, but it can also be used to create more occurences:

.. raw:: html

    <img class='annotation-example'
    src='../_static/images/genbank_annotations/insert_several.png'></img>

This specification can be used both to create new patterns in a region that
contains too few, or to decrease the pattern occurences i a region that contains
too many. Note that with the current algorithm, new occurences of the pattern
will be be preferentially placed towards the center of the selected region.

You can also enforce a sequence (or degenerate sequence) at an exact location
with ``@sequence`` (short for ``@EnforceSequence``):

.. raw:: html

    <img class='annotation-example'
    src='../_static/images/genbank_annotations/enforce_sequence.png'></img>

To enforce several same-length but quite different sequences, use
``@choice`` (short for ``@EnforceChoice``):

.. raw:: html

    <img class='annotation-example'
    src='../_static/images/genbank_annotations/choice.png'></img>

CDS and Codon optimization
--------------------------

To indicate that a region is a CDS and the protein sequence should be conserved
(i.e. only synonymous codon mutations are allowed), use @cds (short for
@EnforceTranslation) on a region whose span is a multiple of 3:

.. raw:: html

    <img class='annotation-example'
    src='../_static/images/genbank_annotations/cds.png'></img>

To codon-optimize a gene use ``~CodonOptimize()``:

.. raw:: html

    <img class='annotation-example'
    src='../_static/images/genbank_annotations/codon_optimize.png'></img>

See `the Codon Tables package webpage <https://github.com/Edinburgh-Genome-Foundry/codon-usage-tables/tree/master/codon_usage_data/tables>`_
for species that can be referred to by name. This includes ``b_subtilis``,
``c_elegans``, ``d_melanogaster``, ``e_coli``, ``g_gallus``, ``h_sapiens``,
``m_musculus``, ``s_cerevisiae``. You can also use a TaxID to refer to a species,
e.g. ``species=1423`` at which case the codon frequencies will be downloaded from
the `Kazusa codon usage database <https://www.kazusa.or.jp/codon/>`_, assuming it not down.

.. caution:: Always use with @cds

   If the CodonOptimize specification is used without a @cds constraint covering
   the same region, then the protein sequence is not guaranteed!

.. caution:: Codon optimization method

    By default, the optimization will replace every codon by the most-frequent codon.
    Use ``method=harmonized_frequencies`` for an optimization method where the
    codons frequency of each amino-acid will reflect the frequencies in the target organism


GC content
----------

Use ``@gc`` to ensure that a given region's GC content is between a
certain range

.. raw:: html

    <img class='annotation-example'
    src='../_static/images/genbank_annotations/gc_range.png'></img>

For large regions, you can use a windowed evaluation, e.g. with the parameter
``window=100`` to ensure that the GC content will remain in the desired range
over every 100bp subsegments of the sequence.

The specification can also be used as an optimization objective, at which case
it is preferable to provide a target rather than a range:

.. raw:: html

    <img class='annotation-example'
    src='../_static/images/genbank_annotations/gc_target.png'></img>

Removing homologies
-------------------

.. raw:: html

    <img class='annotation-example'
    src='../_static/images/genbank_annotations/avoid_non_unique_segments.png'></img>

.. raw:: html

    <img class='annotation-example'
    src='../_static/images/genbank_annotations/avoid_blast.png'></img>

Primer-friendly segments
------------------------

.. raw:: html

    <img class='annotation-example'
    src='../_static/images/genbank_annotations/enforce_melting.png'></img>

.. raw:: html

    <img class='annotation-example'
    src='../_static/images/genbank_annotations/enforce_melting_obj.png'></img>

.. raw:: html

    <img class='annotation-example'
    src='../_static/images/genbank_annotations/allow_primer.png'></img>


Specifications not yet supported as Genbank annotations
--------------------------------------------------------

- AvoidHeterodimerization
- EnforceRegionsCompatibility
