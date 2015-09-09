# coordinates start from 1
SELECT qName,GROUP_CONCAT(DISTINCT feature," ",seqname," ",gff.strand," ",start," ",end)
FROM pseudogene_gff AS gff INNER JOIN ref3_nomis_rep5_psl AS psl
ON gff.strand != '.' AND tName = seqname AND psl.strand = gff.strand AND start < tEnd AND end > tStart+1
GROUP BY qName
