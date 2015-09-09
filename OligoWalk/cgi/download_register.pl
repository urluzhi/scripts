#!/usr/bin/perl  -wT
use CGI qw(:standard);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use Fcntl qw(:flock :seek);
use HTML::Template;
use strict;

##############################################################################
#OligoWalk CGI 								     #
#Copyright 2005 John(Zhi Lu)		zhi_lu@urmc.rochester.edu	     #
#Created Apr. 2, 2008	     #
#Scripts Archives at:			http://rna.urmc.rochester.edu/john   #
##############################################################################

# Get the Date for Entry
my $date =localtime;
#Clear the Env.
$ENV{PATH} = "/bin:/usr/bin:/usr/sbin";
delete @ENV{ 'IFS', 'CDPATH', 'ENV', 'BASH_ENV' };
# open the html template
my $template = HTML::Template->new(filename => '/var/www/html/servers/oligowalk/template.html');
# Set Variables
my $home="http://rna.urmc.rochester.edu";
my $downloadlog = "/var/www/html/servers/oligowalk/download_registration_log.htm";
my $downloadurl = "$home/oligowalk_src.tar.gz";
# Limit the input sequence length
$CGI::DISABLE_UPLOADS = 1;			#not allow upload files 
$CGI::POST_MAX	=1700+40*1024;	 		#The upload sequence size is less than 40K
if (cgi_error()) 	{ &dienice ("The input sequence is too long!"); }


#Get the inputs parameters from web forms, untaint the input 				
my %FORM;
foreach my $p (param()) { 
	if ($p eq 'email') {
		param($p)=~ /([\.\w\-]+)\@([\.\w\-]+)$/;
		$FORM{$p} = "$1\@$2";
	}
	else {
		($FORM{$p}) = param($p) =~/([\w.]+)/ ;
	 
	} 
}


#write log file
&printlog;

# Print Out Output Location Heading 
print "Location: $downloadurl\n\n";



########################################################################################

exit;










#######################
# Subroutines

#error handles
sub dienice {
my($errmsg) = @_;
$template->param(TITLE => "Input Error");
	$template->param(BODY => "<p><span class=\"errormsg\">$errmsg</span></p>");
	# send the obligatory Content-Type and print the template output
	print "Content-Type: text/html\n\n", $template->output; 
    exit;
}


# Log the Entry or Error
sub printlog {
   open (LOG, ">>$downloadlog");
   flock(LOG, LOCK_EX);      # set an exclusive lock 
   seek(LOG, 0, SEEK_END);   # then seek the end of file
   print LOG "[$date] - $FORM{'username'} - $FORM{'institution'} - $FORM{'department'} - $FORM{'email'}<br>\n";
   close (LOG);
}


