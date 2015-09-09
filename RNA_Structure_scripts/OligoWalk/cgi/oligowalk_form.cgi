#!/usr/bin/perl -wT
use CGI;
use CGI qw(:standard);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use HTML::Template;
use strict;

# open the html template
my $template = HTML::Template->new(filename => '/var/www/html/servers/oligowalk/template.html');
#title
$template->param(TITLE => "OligoWalk Web Server for siRNA design");
#filling parts in head
my $head = <<EOF;
<script src="http://rna.urmc.rochester.edu/servers/oligowalk/oligowalk_form.js" type= "text/javascript"> </script>
EOF
$template->param(HEAD => "$head");

#filling body
my $body= <<EOF;

<div id="main">
<div id="head_text">
		    <h1>Welcome to OligoWalk</h1>
				
			<p><strong>OligoWalk</strong> is an online sever calculating thermodynamic features of sense-antisense hybidization. It predicts the free energy changes of oligonucleotides binding to a target RNA. It can be used to <strong>design efficient siRNA</strong> targeting a given mRNA sequence. The source code of OligoWalk for siRNA design can be downloaded from <a href= "http://rna.urmc.rochester.edu/cgi-bin/server_exe/oligowalk/download_form.cgi">here</a>.
<br>The efficient siRNA selection method is described in a published paper (<a href="http://nar.oxfordjournals.org/cgi/content/abstract/36/2/640">link</a>). 
<a href="http://rna.urmc.rochester.edu/servers/oligowalk/help.html#references">More references</a> are listed in the 
<a href="http://rna.urmc.rochester.edu/servers/oligowalk/help.html">Help page</a>. 
</p>
<p>The server has been tested on Firefox 2 and Internet Explorer 7.</p>
</div>

<!-- Testing here: -->
<!--
<form action="/cgi-bin/server_exe/oligowalk/oligowalk.pl" method="post" enctype="multipart/form-data" name="mainform" id="mainform" onsubmit="return" >
-->

<form action="/cgi-bin/server_exe/oligowalk/oligowalk.pl" method="post" enctype="multipart/form-data" name="mainform" id="mainform" onsubmit="return checkinput()" >
<div id=target_text>
				<code> 
					<table>
					<tr >
                        <td colspan="2"> <span class="subheader">Target RNA Sequence</span> </td> 
					</tr>
                    <tr>
                        <td><strong>*Sequence name</strong></td> 
					    <td> 
						   <input name="seqname" type="text" id="seqname" value="example" size="15" maxlength="20" />						</td>
                     </tr>
					 <tr><td><strong>*Primary sequence</strong>
			<a href="javascript:showdiv('seq_id');"><img src="/servers/images/go.gif"   class="my_img" /></a>
					<div class="popup"  id="seq_id">
				<a href="javascript:hidediv('seq_id');"><img src="/servers/images/cancel.jpg" class="icon_img" /></a>
				Each T nucleotide will be converted to U
</div>
					 </td>
				     <td> </td>
					 </tr>
					 <tr>
                        <td>paste the sequence here</td> 
					    <td> 
<textarea name="seqtext" cols="40" rows="300" id="seqtext" >
GCUAUUUUGGUGGAAUUGGUAGACACGAUA
CUCUUAAGAUGUAUUACUUUACAGUAUGAA
GGUUCAAGUCCUUUAAAUAGCACCA
</textarea>				
                     </tr>
                    <tr>
                        <td>or upload the sequence file (file format: fasta)</td> 
					    <td> 
						   <input name="seqfile" type="file" id="seqfile"  />
						</td>
                     </tr>

					 <tr> <td></td><td> </td></tr>
					                        <td><strong>*Email address</strong></td> 
					    <td> 
						  <input name="email" type="text" id="email" size="40" maxlength="50" />
						</td>
                     </tr>
			

	</table>
		
                    </code> 
					<br/>
					</p>	
					<div align="center">
					<input class="button" type="submit" onclick=quicksubmit() /><br/><br/><br/>	
Do not use file uploading, please paste your sequence above for <br/>
 <input class="button"   type="button"  value="Advanced Options"  onclick=next() />
				<br>
				<p>
					</p>
					</div>
					The fields indicated with an asterisk (*) are required. <br><br> 
				</div>
                  <div class=nodisplay style="display:none" id=waiting_text>
				  <p> <strong>Uploading Sequence ... </strong> 
				<br /> (This may take a while)</p>
				  </div>
                    

                  <div class=nodisplay style="display:none" id=option_text>
                   	<code>
					
					 <input class="button_transparent"  name="optiondefault" type="button" id="optiondefault" value="Click to restore the default options of siRNA design"  onclick="option2default()" />
					
                    <table class="mytable">
					<tr >
                        <td colspan="2"><span class="subheader">Oligomer Options</span> </td> 
					</tr>
                      <tr>
                        <td><strong>Oligomer Type</strong></td> 
						<td> 
						<select name="type" size="1" id="type">
                          <option  value="r" selected="selected">RNA</option> 
<!--
						  <option value="d">DNA</option> 
 -->
                        </select>	
					<a href="javascript:showdiv('rnatype');"><img src="/servers/images/go.gif" alt="Only RNA oligo is available now"  class="my_img" /></a>
					<div class="popup"  id="rnatype">
				<a href="javascript:hidediv('rnatype');"><img src="/servers/images/cancel.jpg" class="icon_img" /></a>Only RNA oligo is available now.</div>
						</td>
		
                      </tr>
                      <tr>
                        <td><strong>Oligomer Length</strong></td>
                        <td>
							<input name="length" type="text" id="length" value="19" size="6" maxlength="10"  />						</td>
                      </tr>
                      <tr>
                        <td><strong>Oligomer Concentration</strong></td>
                        <td>
						 <input name="concentration" type="text" id="concentration" value="1" size="6" maxlength="20" />
 					     
                    <select name="unit" size="1" id="unit">
                      <option value="0">M</option>
                      <option value="-3">mM</option>
                      <option value="-6" selected="selected">uM</option>
                      <option value="-9">pM</option>
                      <option value="-12">fM</option>
                    </select>
				<a href="javascript:showdiv('rna_concentration');"><img src="/servers/images/go.gif" alt="Please enter integer number"  class="my_img" /></a>
				<div class="popup"  id="rna_concentration">
				Please enter integer number<a href="javascript:hidediv('rna_concentration');"><img src="/servers/images/cancel.jpg" class="icon_img" /></a>
				</div>
					   </td>
                      </tr>
					  
					  
					  <tr >
                        <td colspan="2">  <span class="subheader">Target Options</span></td> 
												
					  </tr>
					    <tr >
                        <td>  <strong> Prefilter nonfunctional siRNA candidates </strong>
				<a href="javascript:showdiv('prefilter_id');"><img src="/servers/images/go.gif"  class="my_img" /></a>
					<div class="popup"  id='prefilter_id'>
				<a href="javascript:hidediv('prefilter_id');"><img src="/servers/images/cancel.jpg" class="icon_img" /></a>
				Prefilter must be selected if target structure needs to be broken or refolded in 'Binding Mode'.
				</div>

						
						</td> 
						<td> <input name="prefilter" type="checkbox" value="1" checked="checked" />
   
					
						</td>
												
					  </tr>
					    
                    
					  
					  
					  <tr>
                        <td><strong> Binding mode</strong></td>
                        <td>
						  <select name="mode" id="combo_0" onclick="change(this)" onchange="setfoldsize()" >
                          <option value="-1" selected="selected">--Select--</option>
                          <option value="1"> Break local structure</option>
                          <option value="2" > Refold target RNA (default)</option>
                          <option value="3"> Not consider target structure (Fastest)</option>
                          </select>
						</td>
                      </tr>
					  <tr>
                        <td><strong> Folding option</strong></td>
                        <td>
						   <select name="option" id="combo_1" >
                             <option value="-1" selected="selected">--Select--</option>
                             <option value="-2">--Please click on Binding mode's menu again --</option>
                           </select>
						</td>
                      </tr>
					
					  <tr>
                        <td><strong>Folding size</strong></td>
                        <td>
						  <select name="foldsize" id="foldsize"  >
                               <option value="-1" selected="selected">--Select--</option>
                          </select> 
						<a href="javascript:showdiv('foldsize_id');"><img src="/servers/images/go.gif"  class="my_img" /></a>
					<div class="popup"  id='foldsize_id'>
				<a href="javascript:hidediv('foldsize_id');"><img src="/servers/images/cancel.jpg" class="icon_img" /></a>
				The folding size has a limitation of 1000 nucleotides. Please download the software for larger folding size.
				</div>

		     
					     
					      
                     	</td>
                      </tr>
					
					  <tr>
                        <td><strong> Scan region </strong>
						</td>
                        <td>


						   <input name="scanstart" type="text" id="scanstart" value="1" size="5" maxlength="10" >
						   -- <input name="scanstop" type="text" id="scanstop" value="end" size="5" maxlength="10" >
                          
					<a href="javascript:showdiv('scansize');"><img src="/servers/images/go.gif"  class="my_img" /></a>
					<div class="popup"  id="scansize">
				<a href="javascript:hidediv('scansize');"><img src="/servers/images/cancel.jpg" class="icon_img" /></a>The scan region must be less than 1000 nucleotides if siRNA prefilter is not selected. Please download the software if you need longer scan region.</div>


                     	</td>
                      </tr>
					<tr >
                        <td colspan="2"><span class="subheader">Output Options</span> </td> 
					</tr>
					 <tr>
                        <td><strong>Output siRNA candidates selected by SVM</strong></td> 
				       <td> 
						  <input name="siRNA" type="checkbox" value="1" checked="checked" />
				       </td>
                     </tr>
					</table>
                   
                     </code>	
					 </p>           
				    
				
					 <p align="right">
				 	 <input class="button"   type="button"  value="Back"  onclick=goback() />
					 <input class="button" type="submit" />	
					 </p>
					 Use the default if you do not understand the options.<br>
					 The usage of different options are explained in the <a href="http://rna.urmc.rochester.edu/servers/oligowalk/help.html">Help page</a>. <br>

					 </div>



            </form>

		  </div>
		  <div align="left">
		  <!-- Start of StatCounter Code -->
		  <script type="text/javascript">
		  var sc_project=3078565; 
		  var sc_invisible=0; 
		  var sc_partition=27; 
		  var sc_security="bb1aae0d"; 
		  </script>
		  Visitors: 
		  <script type="text/javascript" src="http://www.statcounter.com/counter/counter_xhtml.js"></script><noscript><div class="statcounter"><a class="statcounter" href="http://www.statcounter.com/"><img class="statcounter" src="http://c28.statcounter.com/3078565/0/bb1aae0d/0/" alt="html hit counter" /></a></div></noscript>		  <!-- End of StatCounter Code -->
		  </div>

EOF

# fill in some parameters
$template->param(BODY => "$body");

# send the obligatory Content-Type and print the template output
print "Content-Type: text/html\n\n", $template->output;

