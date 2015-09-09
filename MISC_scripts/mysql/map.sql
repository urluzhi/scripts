SELECT qName,tName,psl.strand,tStart+1,tEnd,GROUP_CONCAT(feature),GROUP_CONCAT(gff.strand),GROUP_CONCAT(start),GROUP_CONCAT(end)
FROM dmel_all_no_analysis_r510_gff AS gff INNER JOIN ref3_nomis_rep5_psl AS psl
ON gff.strand != '.' AND tName = seqname AND psl.strand = gff.strand AND start < tEnd AND end > tStart+1
GROUP BY qName,tName,psl.strand,tStart
