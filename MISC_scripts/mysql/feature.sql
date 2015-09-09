SELECT feature,seqname,gff.strand,start,end,GROUP_CONCAT(qName),GROUP_CONCAT(psl.strand),GROUP_CONCAT(tStart+1),GROUP_CONCAT(tEnd)
FROM pseudogene_gff AS gff INNER JOIN ref3_nomis_rep5_psl AS psl
ON gff.strand != '.' AND tName = seqname AND psl.strand = gff.strand AND start < tEnd AND end > tStart+1
GROUP BY feature,seqname,gff.strand,start
