#!/usr/bin/perl  -wT
use CGI qw(:standard);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use Fcntl qw(:flock :seek);
use HTML::Template;
use strict;

#Clear the Env.
$ENV{PATH} = "/bin:/usr/bin";
delete @ENV{ 'IFS', 'CDPATH', 'ENV', 'BASH_ENV' };


# Set Variables
my $inputhome="/var/www/html/servers/oligowalk/output/";
my $homeurl="http://rna.urmc.rochester.edu";
my $calctimeurl = "$homeurl/servers/oligowalk/help.html#timing";
my $template = HTML::Template->new(filename => '/var/www/html/servers/oligowalk/template.html');

#read input file name
my $input_file;
foreach my $p (param()) { 
	if ($p eq 'file') {
		($input_file)=param($p) =~/([\w.]+)/ ;
		
	}
}
my $outputurl = "$homeurl/cgi-bin/server_exe/oligowalk/oligowalk_out.cgi?file=$input_file";

#test if the input exsits
$input_file =$inputhome.$input_file;
my $content;
unless (-e $input_file) {
$content = <<EOF;
<h1> The calculation is not finished yet.</h1>
<meta http-equiv="refresh" content="10">
<p>This page will be automatically updated in 10 seconds.<br>
The server will also send you an email when your calculation is done.<br>
</p>
<br>Please look at <a href=\"$calctimeurl\" target=\"_blank\">Computation Time Table </a>to estimate the calculation time.<br>
EOF
} 
else {
	open (INPUT, "$input_file");
	while (<INPUT>) {	$content .=$_;	}
	close INPUT;
	unless ($content) {
$content = <<EOF;
<h1> The calculation is not finished yet.</h1>
<meta http-equiv="refresh" content="10">
<p>This page will be automatically updated in 10 seconds.<br>
The server will also send you an email when your calculation is done.<br>
</p>
<br>Please look at <a href=\"$calctimeurl\" target=\"_blank\">Computation Time Table </a>to estimate the calculation time.<br>
EOF

	}
}

#output to template
$template->param(TITLE => "Output of OligoWalk");
$template->param(BODY => "<p>$content</p><br><br>");
print "Content-Type: text/html\n\n", $template->output; 
exit;


print RESULT end_html;
close RESULT;


