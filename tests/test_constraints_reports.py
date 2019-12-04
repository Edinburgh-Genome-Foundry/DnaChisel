import dnachisel as dc
import dnachisel.reports.constraints_reports as cr
import os

# IMPORT THE 10 RECORDS FROM THE genbanks/ FOLDER


def test_constraints_reports():
    genbank_dir = os.path.join("tests", "data", "10_emma_genbanks")
    records = [
        dc.load_record(os.path.join(genbank_dir, filename), name=filename)
        for filename in os.listdir(genbank_dir)
    ]

    # DEFINE THE CONSTRAINTS TO BE CHECKED ON EACH RECORD

    constraints = [
        dc.AvoidPattern("BsaI_site"),
        dc.AvoidPattern("BsmBI_site"),
        dc.AvoidPattern("BbsI_site"),
        dc.AvoidPattern("8x1mer"),
        dc.AvoidPattern("5x3mer"),
        dc.AvoidPattern("9x2mer"),
        dc.AvoidHairpins(stem_size=20, hairpin_window=200),
        dc.EnforceGCContent(mini=0.3, maxi=0.7, window=100),
    ]

    # CREATE A SPREADSHEET AND PLOTS OF THE BREACHES

    dataframe = cr.constraints_breaches_dataframe(constraints, records)
    records = cr.records_from_breaches_dataframe(dataframe, records)
    assert sum([len(r.features) for r in records]) == 157
    pdf_data = cr.breaches_records_to_pdf(records)
    
    assert 70000 < len(pdf_data) < 80000
