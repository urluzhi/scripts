#!/usr/bin/perl  -wT
use CGI qw(:standard);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use Fcntl qw(:flock :seek);
use HTML::Template;
use strict;

##############################################################################
#OligoWalk CGI 								     #
#Copyright 2005 John(Zhi Lu)		zhi_lu@urmc.rochester.edu	     #
#Created Nov. 17, 2005			Last Modified  Dec.2, 2005	     #
#Scripts Archives at:			http://rna.urmc.rochester.edu/john   #
##############################################################################

# Get the Date for Entry
my $date =localtime;

#Clear the Env.
$ENV{PATH} = "/bin:/usr/bin:/usr/sbin";
delete @ENV{ 'IFS', 'CDPATH', 'ENV', 'BASH_ENV' };

# open the html template
my $template = HTML::Template->new(filename => '/var/www/html/servers/oligowalk/template.html');
 


#randomly generate output file names
my $tag =rand(10);		
$tag =~ s/[e\.\+\-]//g;

# Set Variables
my $home="/var/www/html/servers/oligowalk";
my $homeexe="/var/www/cgi-bin/server_exe/oligowalk";
my $homeurl="http://rna.urmc.rochester.edu";
my $reportreal_energy = "$home/output/oligowalk"."$tag"."energy.htm";
my $reportreal_siRNA = "$home/output/oligowalk"."$tag"."siRNA.htm";
my $reportreal_summary = "$home/output/oligowalk"."$tag"."summary.htm";
my $mail_text = "$home/output/oligowalk"."$tag"."mail.txttmp";
my $oligolog = "$home/oligowalk_log.htm";
my $uploadfile ="$home/upload/".$tag.".seq";
my $oligourl ="$homeurl/servers/oligowalk.html";
my $outputurl = "$homeurl/cgi-bin/server_exe/oligowalk/oligowalk_out.cgi?file=oligowalk".$tag."siRNA.htm";
my $outputurl_energy = "$homeurl/cgi-bin/server_exe/oligowalk/oligowalk_out.cgi?file=oligowalk".$tag."energy.htm";
my $outputurl_summary = "$homeurl/cgi-bin/server_exe/oligowalk/oligowalk_out.cgi?file=oligowalk".$tag."summary.htm";
my $calctimeurl = "$homeurl/servers/oligowalk/help.html#timing";
my $job = "/var/www/cgi-bin/server_exe/oligowalk/jobs4sge/oligo".$tag.".q";
my $recipient = 'David_Mathews@urmc.rochester.edu';


# Set Your Options:
my $redirection = 0;       # 1 = Yes; 0 = No
my $remote_mail = 0;       # 1 = Yes; 0 = No; not sending by perl because of permission problem
my $mail_by_cron = 2;      # 1 = html format by sendmail; 2= text format by mailsend(smtp);  0 = No; sending by cron job
 

# If you answered 1 to  $remote_mail you will need to fill out 
# these variables below:
my $mailprog = 'sendmail';

# Done
##############################################################################
#check if queue is full
&checkqueue;

# Limit the input sequence length
$CGI::DISABLE_UPLOADS = 0;			#allow upload files 
$CGI::POST_MAX	=1700+40*1024;	 		#The upload sequence size is less than 40K
if (cgi_error()) 	{ &dienice ("The input sequence is too long!"); }


#Get the inputs parameters from web forms, untaint the input 				
my %FORM;
foreach my $p (param()) { 
	if ($p eq 'email') {
		param($p)=~ /([\.\w\-]+)\@([\.\w\-]+)$/;
		$FORM{$p} = "$1\@$2";
	}
	elsif($p eq 'seqtext') {
		($FORM{$p}) = param($p) =~/([\cM\s\w.]+)/ ;

	}
	elsif ($p eq 'seqfile'){
		 $FORM{$p}=param($p);
	}
	else {
		($FORM{$p}) = param($p) =~/([\w.]+)/ ;
	 
	} 
#	unless ( $FORM{$p} ) {	&dienice ("Unreadable input, Please check your input format.<br>\n");}
}
$FORM{'unit'}= -$FORM{'unit'};
#check the validation of email address
my $valid = 0;
if ($FORM{'email'} =~ /^\w+[\.\w\-]+\@[\.\w\-]+\.[\w\-]+$/) {$valid = 1;}
&no_email unless ($valid);

#check the options
if ($FORM{'option'} == 2 && $FORM{'mode'} != 2)		{ &wrong_option; }
if ($FORM{'foldsize'} ne 'all' && $FORM{'mode'} ==3)  { &wrong_foldsize; }
if ($FORM{'option'} == 1 && $FORM{'mode'} == 3)		{ &dienice("Option 1 cannot be chosen when Mode 3\n"); }

#read the sequence  and count the number of nucleotides
my $rna;
if (-e $FORM{'seqfile'}) {#read from uploaded file
	&dienice ("The uploaded sequence file is not a text file!\n")	unless(-T $FORM{'seqfile'});
	my $filehandle = $FORM{'seqfile'};
	my @filedata =<$filehandle>;	#this will store the file content
	$rna = extract_sequence_from_fasta_data(@filedata);
}
else { # read from pasted text
	$rna=$FORM{'seqtext'};
}
#check the rna sequence:
if ($rna eq '')		{&dienice("Sequence file has wrong format\n"); }
$rna=~ s/\cM//g;
$rna=~ s/\s//g;
$rna=~ s/[0-9]//g;
$rna=~ tr/acgtuT/ACGUUU/;  # convert each T to U
my $checkrna=$rna;
$checkrna=~ s/[AGTCUagtcuxXN]//g;
&dienice("Sequence file or pasted sequence has wrong format, e.g. characters other than AGTCUXNagtcux were found. Please note that the readable file format is fasta.\n")	unless($checkrna eq '');  

my @rna= split('',$rna);
my $num=@rna;		#number of sequence nucleotides

#read the interger number of concentration
$FORM{'concentration'}= int($FORM{'concentration'});

#check input errors
my $checktarget = $num ;
my $checkfold   = $FORM{'foldsize'} ;
if ($checkfold eq 'all') 	{ $checkfold =$checktarget; }
if( $checktarget  > 10000 ) {
	&dienice ("Target sequence should be less than 10000 nt for the server.\n Please go to the homepage and download the source code if you want large calculation.\n");
} 
if( $checkfold  > 1500 ) {
	&dienice ("Folded region should not be longer than 1500 nt for the server\n Please go to the homepage and download the source code if you want larget calculation.\n");
} 
if( $FORM{'length'}  >= $checkfold  ) {
	&dienice ("Oligomer seems to be too long or folding size is too small\n");
}

# Print Out Output Location Heading 
if ($redirection eq '1') {
   print "Location: $oligourl\n\n";
}
else { 
	&no_redirection;	#print the heading of output
}

# email the output link to users from perl; it is not used
if ($remote_mail eq '1' && $FORM{'email'}) {
   open (MAIL, "|$mailprog -t") || &dienice("Can't open $mailprog!\n");

   print MAIL "To: $FORM{'email'}\n";
   print MAIL "Cc: $recipient\n";
   print MAIL "From: $recipient\n";
   print MAIL "Subject: OligoWalk Calculation Complete\n\n";
   print MAIL "Thank you for using OligoWalk.\n\n";
   print MAIL "Please get your result from:\n $outputurl_summary \n";
   print MAIL "------------------------------------------------------\n\n\n";

   close (MAIL);
}
#alternative email method: writing simple text for mailing
if ($mail_by_cron eq '1' && $FORM{'email'}) {
   open (MAIL, "> $mail_text") || &dienice("Can't open $mail_text !\n");
   print MAIL "To: $FORM{'email'}\n";
   print MAIL "Cc: $recipient\n";
   print MAIL "From: $recipient\n";
   print MAIL "Subject: OligoWalk Calculation Complete\n";
   print MAIL <<EOF;
Content-type: text/html

<HTML> 
<HEAD> <TITLE>Thank you</TITLE></HEAD>
<BODY> 
<h2>Thank you for using OligoWalk web server</h2>
<br><HR>
The OligoWalk calculation is finished.<br>
<a href="$outputurl_summary">Click here to see the report.</a><br><br><br>
<HR>
Please email questions or errors to $recipient with reference number: $tag . <br>
Thanks. <br> <br>
</BODY> 
</HTML> 

EOF
close (MAIL);

}
elsif ($mail_by_cron eq '2' && $FORM{'email'}) {
   open (MAIL, "> $mail_text") || &dienice("Can't open $mail_text !\n");
   print MAIL "smtp.urmc.rochester.edu\n";
   print MAIL "urmc.rochester.edu\n";
   print MAIL "$recipient\n";
   print MAIL "$FORM{'email'}\n";
   print MAIL <<EOF;
Thank you for using OligoWalk web server.
-----------------------------------------
The OligoWalk calculation is finished.
The link of the result is
$outputurl_summary

-----------------------------------------
Please email questions or errors to $recipient with reference number: $tag . 
Thanks.

EOF
close (MAIL);

}


#upload the file to be read by oligowalk.o 
open (UPLOADFILE, ">$uploadfile") or  &dienice("Can't open .seq file: $!\n");
flock(UPLOADFILE, LOCK_EX);      # set an exclusive lock 
seek(UPLOADFILE, 0, SEEK_SET);   # then seek the begining of file
print UPLOADFILE  	";\n";
$FORM{'seqname'} =~ s/\n//g;
print UPLOADFILE	"$FORM{'seqname'}\n";
print UPLOADFILE	"$rna";
print UPLOADFILE	"1\n"; 
close UPLOADFILE; 

#######################################################################################
#initiated the value which is readable by oligowalk.o 
if (!$FORM{'prefilter'})		{ $FORM{'prefilter'}=0; }
if ($FORM{'scanstop'} eq 'end')		{ $FORM{'scanstop'}=0; }
if ($FORM{'foldsize'} eq 'all')		{ $FORM{'foldsize'}=0; }
else					{ $FORM{'foldsize'}=$FORM{'foldsize'} - $FORM{'length'}; }
my $mainstring="-type  $FORM{'type'} -seq  $uploadfile -o nofile -m $FORM{'mode'} -st 1 -en 0 -M  $FORM{'scanstart'} -N  $FORM{'scanstop'} -l  $FORM{'length'} -s  $FORM{'option'} -co  $FORM{'concentration'} -unit  $FORM{'unit'} -fi  $FORM{'prefilter'} -fold  $FORM{'foldsize'} -score ";


#writing job srcipt for sge
&sending_job ($mainstring);

########################################################################################

exit;










#######################
# Subroutines

#error handles
sub dienice {
my($errmsg) = @_;
$template->param(TITLE => "Input Error");
	$template->param(BODY => "<p><span class=\"errormsg\">$errmsg</span></p><p><br></p><p>Return to <a href=\"$oligourl\">OligoWalk</a>. </p><p><br></p>");
	# send the obligatory Content-Type and print the template output
	print "Content-Type: text/html\n\n", $template->output; 
    exit;
}

sub no_email {
	# fill in some parameters
	$template->param(TITLE => "Email Address Error");
	$template->param(BODY => "<p><h1>Please input the correct email address</h1></p><p><br></p><p>Return to <a href=\"$oligourl\">OligoWalk</a>. </p><p><br></p>");
	# send the obligatory Content-Type and print the template output
	print "Content-Type: text/html\n\n", $template->output; 
	exit;
}


sub wrong_option {
	print header;
	print start_html("Sorry: Wrong option");
	print "<h1>Sorry: Wrong option</h1>\n"; 
	print "<p>Option 2 (partition function) can only be chosen when Mode 2 (refold target) ";
	print "is used.<p>\n";
	print "Return to the <a href=\"$oligourl\">OligoWalk</a>.";
	print end_html;
	exit;
}

sub wrong_foldsize {
	print header;
	print start_html("Sorry: Folding size cannot be redefined");
	print "<h1>Sorry: Folding size cannot be redefined</h1>\n"; 
	print "<p>Folding size can only be redefined when Option 2 (partition function) ";
	print "is used.<p>\n";
	print "Return to the <a href=\"$oligourl\">OligoWalk</a>.";
	print end_html;
	exit;
}


# Log the Entry or Error
sub log {
   my($log_type) = @_;
   open (LOG, ">>$oligolog");
   flock(LOG, LOCK_EX);      # set an exclusive lock 
   seek(LOG, 0, SEEK_END);   # then seek the end of file
   if ($log_type == 0) {
      print LOG "[$date] - $FORM{'email'} - $tag - Job submitted.<br>\n";
   }
   else{
      print LOG "[$date] - $FORM{'email'} - $tag - Failed to submit the job. <br>\n ";
   }
   close (LOG);
}

# Print out the heading of output
sub no_redirection {


	my $content= <<EOF;
<h1> Thank you for using OligoWalk.</h1> <br>
Job submitted on 
$date <br>

EOF
	if ($FORM{'siRNA'}) {
		$content .= <<EOF;
<b><a href=\"$outputurl\" target=\"_blank\">Click here </a></b> to get the siRNA candidates <br> 
<b><a href=\"$outputurl_energy\" target=\"_blank\">Click here </a></b> to get the free energy table 
</p>
<br><b>Options:</b><hr>

EOF
	}
	else {
		$content .= <<EOF;
<b><a href=\"$outputurl_energy\" target=\"_blank\">Click here </a></b> to get the free energy table 
</p>
<br><b>Options:</b><hr>

EOF
	}

	if ($FORM{'prefilter'}) {
		$content .= "prefilter is used<br>\n";
	}
	else {
		$content .= "prefilter is not used<br>\n";
	}
	if ($FORM{'option'} == 0) {
		$content .= "option = optimal structure prediction for target<br>\n";
	}
	elsif ($FORM{'option'} == 3) {
		$content .= "option = suboptimal structure prediction for target<br>\n";
	}
	elsif ($FORM{'option'} == 2) {
		$content .= "option = partition function calculation for target<br>\n";
	}
	if ($FORM{'mode'} == 1) {
		$content .= "mode = breaking local structure of target<br>\n";
	}
	elsif ($FORM{'mode'} == 2) {
		$content .= "mode = refold target RNA<br>\n";
	}
	elsif ($FORM{'mode'} == 3) {
		$content .= "mode = not consider target structure<br>\n";
	}
	foreach  my $name (keys %FORM) {
      if ($name eq 'unit')	{ next; }
      if ($name eq 'foldsize' and $FORM{'mode'}==3)	{ next; }
      if ($name eq 'prefilter') { next; }
      if ($name eq 'option') { next; }
      if ($name eq 'mode') { next; }
      if ($name eq 'siRNA') { next; }
      if ($name eq 'seqtext')   { next; }
      if ($name eq 'concentration') {
		  $content .= "$name = $FORM{'concentration'} e $FORM{'unit'} (M)<br>\n"; 
      }
      else	{ $content .= "$name = $FORM{$name}<br>\n";}
   }
	if ($FORM{'siRNA'}) { $content .="siRNA = Yes<br>\n";	}
   #print out the sequence with multiple lines
   $content.= "<br><b>Input sequence</b> ( $num bases):<br>\nEach base T has been converted to base U<br><HR>\n";
   for (my $pos=0; $pos < length($rna); $pos += 20) {
   	$content .= substr($rna, $pos, 20)."<br>\n";
   }
   $content.= "<br>";
   $content .= <<EOF;
<h2> The server will send you an email with this link when your calculation is done.</h2> <br> 

EOF

	$template->param(TITLE => "Submitting Job");
	$template->param(BODY => "<p>$content</p><br><br>");
	print "Content-Type: text/html\n\n", $template->output;
	#output to a summary text file
	open (FILE, ">$reportreal_summary") or &dienice ("cannot open $reportreal_summary"); 
	print FILE $content;
	close (FILE);
}

sub sending_job {
	my ($input_options)=@_;
	my $tmp= <<EOF;
#!/bin/bash
#\$ -S /bin/bash
#submit jobs from current directory
#\$ -cwd
#wall time
#\$ -l h_rt=1000:00:00
#standard output file, error also output to the same file
#\$ -o $job.out -j y

DATE1=`date`
export DATAPATH=/var/www/cgi-bin/server_exe/NN_data/
source /var/www/cgi-bin/server_exe/sge_root/default/common/settings.sh
echo  "Start at \$DATE1"

oligowalk.exe $input_options > $reportreal_energy

DATE2=`date`
echo  "Finish at \$DATE2"
echo "<p>Start at \$DATE1 <br />\nFinish at \$DATE2 </p>" >> $reportreal_energy 
echo "<p>Start at \$DATE1 <br />\nFinish at \$DATE2 </p>" >> $reportreal_summary 
EOF

	if ($FORM{'siRNA'}) {
		$tmp .= "cat $reportreal_energy | /var/www/cgi-bin/server_exe/oligowalk/oligowalk_post.pl $reportreal_siRNA \n";
	}
	$tmp .= "rename txttmp txt $mail_text \n ";
	open (SGE,">$job");
	print SGE $tmp,"\n";
	close SGE;
	my $main_return = system("source /var/www/cgi-bin/server_exe/sge_root/default/common/settings.sh; qsub $job");
   	&log($main_return);

}

sub checkqueue {
	open (QUEUE, 'source /var/www/cgi-bin/server_exe/sge_root/default/common/settings.sh; qstat  |') ;
	my @lines= <QUEUE>;
	close QUEUE;
	if (scalar @lines > 12) { # This means the maximum jobs in the queue is N-2 
		dienice("Server busy,Please resubmit later");
	}
}


# A subroutine to extract FASTA sequence data from an array
# adapted from " BeginPerlBioinfo.pm"
sub extract_sequence_from_fasta_data {
    my(@fasta_file_data) = @_;

	# Declare and initialize variables
    my $sequence = '';

    foreach my $line (@fasta_file_data) {

        # discard blank line
        if ($line =~ /^\s*$/) {
            next;

        # discard comment line
        } elsif($line =~ /^\s*#/) {
            next;

        # discard fasta header line
        } elsif($line =~ /^>/) {
            next;

        # keep line, add to sequence string
        } else {
            $sequence .= $line;
        }
    }

    # remove non-sequence data (in this case, whitespace) from $sequence string
    $sequence =~ s/\s//g;

    return $sequence;
}

