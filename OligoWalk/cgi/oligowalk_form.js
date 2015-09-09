// JavaScript Document
<!--
//count the bases
var numofbases;
var seq=document.mainform.seqtext.value;
//remove white space in textarea
seq = seq.replace(/\s/g,""); 
numofbases=seq.length;

//control options and modes 	
function change(currentbox) {	
//  mode box
//  mode_1 = new Option("(1) Break local structure", "1");
//	mode_2 = new Option("(2) Refold target RNA (default)", "2");
//	mode_3 = new Option("(3) Not consider target structure (Fastest)", "3");
//  option box

	option_0 = new Option("--Select--", "-1");
	option_1 = new Option(" Only consider the Optimal Structure", "0");
	option_2 = new Option(" Consider Suboptimal Structures", "3");
	option_3 = new Option(" Consider all possible structures using partition function", "2");

	son = document.mainform.option;
	// I reset the first option
	for (m=son.options.length-1;m>0;m--) son.options[m]=null;
	//decide the options number 
	m=-2;
	if (currentbox.value == "-1" ){		
		son.options[0]=new Option("--Select--","-1");
	}
	else if(currentbox.value == "1" ){		    
		m=1;
	}
	else if ( currentbox.value == "2" ){		
		m=2;
	}
	else if ( currentbox.value == "3" ){
		son.options[0]=new Option("No Options for this binding mode","99");
	}
	// filling the "son" combo (if exists)
	for (i=m+1;i>=0;i--){
    	eval("son.options["+i+"]=new Option(option_"+i+".text, option_"+i+".value)");
	}
	if (m==1 || m==2) {	son.options[2].selected =true;}
	else   {son.options[0].selected =true;}
		
}

//default settings for option and mode
function option2default() {
	option_0 = new Option("--Select--", "-1");
	option_1 = new Option(" Only consider the Optimal Structure", "0");
	option_2 = new Option(" Consider Suboptimal Structures", "3");
	option_3 = new Option(" Consider all possible structures using partition function", "2");
	son = document.mainform.option;
	document.mainform.mode.options[2].selected=true;
	for (i=3;i>=0;i--)    	eval("son.options["+i+"]=new Option(option_"+i+".text, option_"+i+".value)");
	son.options[3].selected=true;
	document.mainform.type.value='r';
	document.mainform.length.value=19;
	document.mainform.prefilter.checked=true;
	document.mainform.prefilter.checked=true;
	document.mainform.siRNA.checked=true;
	setfoldsize();	
}

//check if the uploading field is blank for Advanced Options
function checkupload () {
	if ( document.mainform.seqfile.value != "") {
		alert ("Please clear your file uploading field if you want Advanced Options.\n Just paste the primary sequence to the text area.");
		return false;
	}
	else {
		return true;
	}
}



//check the seqtext format and length
function checkseqtext() {
	var x=document.mainform.seqtext.value;
	countbases();
	if ( numofbases > 10000) {
		alert(" The maximum sequence length is 10000.Please truncate the sequence. ");
		return false;
	}

	var bases=/^[aguctAGUCT \r\n\t]+$/;
	if (x.match(bases) != null){
		testresult=true;
	}
	else{
		alert("Sequence format Error:\n Invalid base(s) found, Please input the correct base:\n a g u c t A G U C T");
		testresult=false;
	}
	return (testresult);
}

//check the options
function checkoptions() {
	var testresult;
	
	if (document.mainform.mode.value <0 ){ 
		alert("Please select a valid Binding mode!");
		testresult=false;
	}
	else if (document.mainform.option.value <0 ){ 
		alert("Please select a valid Folding option!");
		testresult=false;
	}
	else if (document.mainform.mode.value ==3 ) {
		document.mainform.foldsize.value = "all";	
		testresult=true;
	}
	else {
		testresult=true;
	}
	
	if (document.mainform.mode.value != 3 && document.mainform.prefilter.checked != true ) {
		alert("Prefilter must be selected for this binding mode.\nYou can run the program without prefilter by downloading the OligoWalk program.");
		testresult=false;
	}

	return (testresult);
}

//check oligomer
function checkoligo() {
	var testresult;	
	var numbers=/(^\d+$)|(^\d+\.\d+$)/;
	if ( !numbers.test(document.mainform.concentration.value) || !numbers.test(document.mainform.length.value) ) {
		alert ("Oligomer Error:\n Please input positive numbers for length and concentration values.");						
		testresult=false;
	}
	else if ( document.mainform.concentration.value > 100000 || document.mainform.concentration.value < 1){ 
		alert ("Oligomer Error:\n The value of concentration is out of rang\n choose larger unit or smaller number");
		testresult=false;
	}
	else if (document.mainform.length.value > 500 || document.mainform.length.value < 1 )  {
		alert("Oligomer Error:\n The length of the oligo is out of range (1-500)");
		testresult=false;
	}
	else {
		testresult=true;
	}
	return (testresult);
}

//check target
function checktarget() {
	var testresult;	
	var numbers=/(^\d+$)|(^\d+\.\d+$)/;
	var oligolength=document.mainform.length.value;
	var scanstart =document.mainform.scanstart.value;
	var scanstop =document.mainform.scanstop.value;
	var foldsize =document.mainform.foldsize.value;
	if (scanstop == 'end') {	
		scanstop=numofbases-oligolength+1;
	}
	if (foldsize == 'all') {	
		foldsize=numofbases;
	}

	if ( !numbers.test(foldsize) || !numbers.test(scanstart) || !numbers.test(scanstop) ) {
		alert ("Target Error:\n Please input positive numbers for scan region.");						
		testresult=false;
		return false;
	}
	foldsize=parseInt(foldsize);
	if (foldsize == 'all' ){	
		foldsize=numofbases;
	}
	oligolength= parseInt(oligolength);
	scanstop= parseInt(scanstop);
	scanstart= parseInt(scanstart);
	var scanlength;
	scanlength = scanstop - scanstart +1;
	if ( foldsize <0) {
		alert(" Please select folding size again.");
		testresult=false;
	}
	else if (  foldsize > 1000	&& document.mainform.mode.value !=3 ) {
		alert(" The maximum folding size is 1000.Please download the software if you need larger folding size. We recommend 800. ");
		testresult=false;
	}
	else if ( numofbases > 10000) {
		alert(" The maximum sequence length is 10000.Please truncate the sequence. ");
		testresult=false;
	}

	else if ( scanstart <0 || scanstop <0) {
		alert(" Please select scanning region again.");
		testresult=false;
	}
	else if ( foldsize > numofbases || foldsize <= oligolength){ 
		alert ("Taget Error:\n The folding size or oligomer length is out of range.\n It must larger than the length of oligomer and less than the whole length of target.");
		testresult=false;
	}
	else if ( scanlength <0 ) {
		alert ("Start site should be less than the end.");
		testresult=false;
	}
	else if ( document.mainform.prefilter.checked != true && scanlength > 1000 ) {
		alert ("The scan region must be less than 1000 nucleotides if prefilter is not chosen.Please download the software for more scanning.");
		testresult=false;
	}

	else if ( scanstop > numofbases - oligolength +1 ) {
		alert ("Target Error:\nThe end of scan region should be less than or equall to \n target length - oligo length + 1.");
		testresult=false;
	}
	else {
		testresult=true;
	}
	return (testresult);
}

//validation of email address
function checkemail(){
	var testresults;
	var str=document.mainform.email.value;
	var filter=/^([\w-]+(?:\.[\w-]+)*)@((?:[\w-]+\.)*\w[\w-]{0,66})\.([a-z]{2,6}(?:\.[a-z]{2})?)$/i;
	if (filter.test(str)) {
		testresults=true;
	}
	else{
		alert("Email Address Error:\nPlease input a valid email address!");
		testresults=false	;
	}
	return (testresults);
}

//Block multiple form submission script
function checksubmit(){
	document.mainform.submit();	
	checksubmit=blocksubmit;
	return false;
}
function blocksubmit(){
	var formerrormsg="You\'ve attempted to submit the form multiple times.\n Please reload page if you need to resubmit form.";
	if (typeof formerrormsg!="undefined") 		alert(formerrormsg);
	return false;
}

//count the bases number
function countbases() {
	var seq=document.mainform.seqtext.value;
	//remove white space in textarea
	seq = seq.replace(/\s/g,""); 
	numofbases=seq.length;
}

//show the foldsize option
function listfoldsize() {
	var foldsize=document.mainform.foldsize;

	for (m=foldsize.options.length-1;m>0;m--) foldsize.options[m]=null;
	if (document.mainform.mode.value ==3) {
		foldsize.options[0]=new Option("No folding","all");
		foldsize.options[0].selected = true;

	}
	else {
		// filling the the foldsize menu (if exists)
		m=0;
		var foldsize_options;
		if (numofbases > 1000) {
			foldsize_options =1001;
		}
		else {
			foldsize_options =numofbases;
		}
		for (i=1;i<foldsize_options/100;i++){
			m++;
			eval("foldsize.options["+m+"]=new Option("+m*100+", "+m*100+")");
		}	
		foldsize.options[0]=new Option("--Select--","-1");
		if (numofbases <= 1000) {
			foldsize.options[m+1]=new Option("all","all");
		}
		if (numofbases > 800) {
			foldsize.options[8].selected = true;
					}
		else {
			foldsize.options[m+1].selected = true;
		}	 
	}
}
//set the foldsize option
function setfoldsize() {
	var foldsize=document.mainform.foldsize;
	if (document.mainform.mode.value ==3) {
		foldsize.options[0]=new Option("No folding","all");
		foldsize.options[0].selected = true;

	}
	else {
		foldsize.options[0]=new Option("--Select--","-1");
		if (numofbases > 800) {
			foldsize.options[8].selected = true;
					}
		else {
			foldsize.options[foldsize.options.length-1].selected = true;
		}	 
 
	}
}

//show the scan region's options
function listscansize() {
	var scanstart=document.mainform.scanstart;
	var scanstop=document.mainform.scanstop;
	var length=document.mainform.length.value;
	//reset all options except the first one
	for (m=scanstart.options.length-1;m>0;m--) scanstart.options[m]=null;
	for (m=scanstop.options.length-1;m>0;m--) scanstop.options[m]=null;

	// filling the menu 
	for (i=numofbases-length;i>0;i--)    	scanstart.options[i]=new Option(i, i);
	for (i=numofbases-length;i>0;i--)    	eval("scanstop.options["+i+"]=new Option("+i+", "+i+")");
	for (i=numofbases-length+1;i>numofbases-length;i--)    	scanstop.options[i]=new Option("end",i);
	scanstart.options[0]=new Option("--Select--","-1");
	scanstop.options[0]=new Option("--Select--","-1");
	scanstart.options[1].selected=true;
	scanstop.options[numofbases-length+1].selected=true;
}

//reset foldsize when option changed 
function resetfoldsize() {
	var foldsize=document.mainform.foldsize;
	if (foldsize.value != "0") {
		//reset all options except the first one
		for (m=foldsize.options.length-1;m>0;m--) foldsize.options[m]=null;
		foldsize.options[0]=new Option("--Select--","-1");
	}
}

//reset foldsize 
function resetsize() {
	resetfoldsize();
}

//check the inputs and submit the job only once
function checkinput() {
	return (checkemail() &&  checkoptions() && checkoligo() && checktarget() && checksubmit() );
}

/* This tells the browser to show just one id while hiding all the others. */

function next() {
	hidediv('target_text');
//	showdiv('waiting_text');
	setTimeout('listing()',10);
	option2default();
}

function quicksubmit() {
	if (checkseqtext()) { 	listfoldsize();	}
	option2default();
}

function listing() {
	if (checkseqtext() && checkemail() && checkupload() ) {
	listfoldsize();
	hidediv('waiting_text');
	hidediv('head_text');
	showdiv('option_text');
	}
	else {
	hidediv('waiting_text');
	showdiv('head_text');
	showdiv('target_text');
	}	
}

function goback() {
		hidediv('option_text');
		showdiv('head_text');
		showdiv('target_text');
}

function hidediv(id) {
	//safe function to hide an element with a specified id
	if (document.getElementById) { // DOM3 = IE5, NS6
		document.getElementById(id).style.display = 'none';
	}
	else {
		if (document.layers) { // Netscape 4
			document.id.display = 'none';
		}
		else { // IE 4
			document.all.id.style.display = 'none';
		}
	}
}

function showdiv(id) {
	//safe function to show an element with a specified id
		  
	if (document.getElementById) { // DOM3 = IE5, NS6
		document.getElementById(id).style.display = 'block';
	}
	else {
		if (document.layers) { // Netscape 4
			document.id.display = 'block';
		}
		else { // IE 4
			document.all.id.style.display = 'block';
		}
	}
}


//-->	


