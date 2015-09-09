CREATE TABLE `ref3_nomis_rep5_psl` (
# `id` int unsigned auto_increment,
  `matches` tinyint(2) unsigned NOT NULL default '0',
  `misMatches` tinyint(2) unsigned NOT NULL default '0',
  `repMatches` tinyint(2) unsigned NOT NULL default '0',
  `nCount` tinyint(2) unsigned NOT NULL default '0',
  `qNumInsert` tinyint(2) unsigned NOT NULL default '0',
  `qBaseInsert` tinyint(2) unsigned NOT NULL default '0',
  `tNumInsert` tinyint(2) unsigned NOT NULL default '0',
  `tBaseInsert` tinyint(2) unsigned NOT NULL default '0',
  `strand` char(1) NOT NULL default '',
  `qName` varchar(255) NOT NULL default '',
  `qSize` tinyint(2) unsigned NOT NULL default '0',
  `qStart` tinyint(2) unsigned NOT NULL default '0',
  `qEnd` tinyint(2) unsigned NOT NULL default '0',
  `tName` varchar(255) NOT NULL default '',
  `tSize` int(10) unsigned NOT NULL default '0',
  `tStart` int(10) unsigned NOT NULL default '0',
  `tEnd` int(10) unsigned NOT NULL default '0',
  `blockcount` tinyint(2) unsigned NOT NULL default '0',
  `blockSizes` varchar(50) NOT NULL default '',
  `qStarts` varchar(50) NOT NULL default '',
  `tStarts` varchar(50) NOT NULL default '',
   
   KEY `qName` (`qname`(7),`tStart`)

) ;


