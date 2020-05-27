rule qpcr_temperature_effect:
    input:
        cq = "data/qpcr/190514/190514_bl6_coldIBat.txt",
        sample_info = "sample_info/sampleInfo.csv"
    output:
        "res/qpcr/qpcr_temperature_effect.html"
    script:
        "qpcr_temperature_effect.Rmd"

rule qpcr_dtu_validation:
    input:
        cq1 = "data/qpcr/200514/200514_bl6_dtu_1.txt",
        cq2 = "data/qpcr/200514/200514_bl6_dtu_2.txt",
        sample_info = "sample_info/sampleInfo.csv",
        signif = "res/comparisons/comparisons_dtu_significant.csv"
    output:
        "res/qpcr/qpcr_dtu_validation.html"
    script:
        "qpcr_dtu_validation.Rmd"
