CREATE TABLE `pseudogene_gff` (
  `seqname` varchar(255) NOT NULL default '',
  `source` varchar(255) NOT NULL default '',
  `feature` varchar(255) NOT NULL default '',
  `start` int(10) unsigned NOT NULL default '0',
  `end` int(10) unsigned NOT NULL default '0',
  `score` smallint(4) unsigned NOT NULL default '0',
  `strand` char(1) NOT NULL default '',
  `frame` char(1) NOT NULL default '',
  `group` varchar(255) NOT NULL default '',

   KEY `seqname` (`seqname`(7),`start`)

) ;


