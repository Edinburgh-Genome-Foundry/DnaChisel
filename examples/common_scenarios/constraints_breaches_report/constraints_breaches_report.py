import dnachisel as dc
import dnachisel.reports.constraints_reports as cr
import os

# IMPORT THE 10 RECORDS FROM THE genbanks/ FOLDER

records = [
    dc.load_record(os.path.join("genbanks", filename), name=filename)
    for filename in os.listdir("genbanks")
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
dataframe.to_excel("breaches.xlsx")
records = cr.records_from_breaches_dataframe(dataframe, records)
cr.breaches_records_to_pdf(records, "breaches_plots.pdf")

print("Done! Check breaches.xlsx and breaches_plots.pdf for results.")