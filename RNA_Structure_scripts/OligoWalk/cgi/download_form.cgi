#!/usr/bin/perl -wT
use CGI;
use CGI qw(:standard);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use HTML::Template;
use strict;

# open the html template
my $template = HTML::Template->new(filename => '/var/www/html/servers/oligowalk/template.html');
#title
$template->param(TITLE => "OligoWalk Registration Form");
#filling parts in head
my $head = <<EOF;
<script src="http://rna.urmc.rochester.edu/servers/oligowalk/oligowalk_form.js" type= "text/javascript"> </script>
EOF
$template->param(HEAD => "$head");

#filling body
my $body= <<EOF;

<h1>OligoWalk Registration Form for Downloading</h1>
The information is only used as a record of our server. It will not be used for any other purpose. Thank you.				
<form action="/cgi-bin/server_exe/oligowalk/download_register.pl" method="post" enctype="multipart/form-data" name="mainform" id="mainform" onsubmit="return checkemail()" >
<code> 
	<table>

	<tr>
    <td><strong>Name:</strong></td> 
	<td> <input name="username" type="text"  size="15" maxlength="50" /></td>
    </tr>
	
	<tr>
    <td><strong>Institution/University:</strong></td> 
	<td> <input name="institution" type="text"  size="15" maxlength="50" /></td>
    </tr>
	
	<tr>
    <td><strong>Department:</strong></td> 
	<td> <input name="department" type="text"  size="15" maxlength="50" /></td>
    </tr>
	
	<tr>
    <td><strong>E-mail:</strong></td> 
	<td> <input name="email" id="email" type="text"  size="15" maxlength="50" /></td>
    </tr>


	</table>
	<p> 	<input class="button" value="Download Now" type="submit" />	</p>
</code> 
EOF

# fill in some parameters
$template->param(BODY => "$body");

# send the obligatory Content-Type and print the template output
print "Content-Type: text/html\n\n", $template->output;

