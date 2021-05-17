version 1.0
# pzw
# 20210517
# Lung Cancer Tumor-only mode(FFPE)

import "../Task/QC/QC.wdl" as qc
import "../Task/Mapping/Mapping.wdl" as mapping

workflow Lung {
    input {
        String sample
        File rawRead1
        File rawRead2
        Int threads
    }

    call qc.Fastp as QC {
        input:
            sample = sample,
            threads = threads,
            rawRead1 = rawRead1,
            rawRead2 = rawRead2
    }

    call qc.FastpFormat as FastpReport {
    	input:
    		sample = sample,
    		jsonReport = QC.jsonReport
    }

    call mapping.Bwa as Mapping {
        input:
            sample = sample,
            threads = threads,
            cleanRead1 = QC.cleanRead1,
            cleanRead2 = QC.cleanRead2
    }
}
