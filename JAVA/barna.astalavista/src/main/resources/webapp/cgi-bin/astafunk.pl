#!/usr/bin/perl
print "Content-type: text/html\n\n";
print "<h1>Hello World</h1>\n";

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);

#use File::Temp qw/tempfile/;
use File::Temp qw/tempdir/;

my $_debug = 1;

my $APACHE_ROOT = $ENV{'APACHE_ROOT'};
my $HOME = "~gmaster";
### MICHA ###
# changed by 22.01.08, $HOME points via apache to the old home.
# $ENV{'HOME'};
my $FRAME  = "$APACHE_ROOT/htdocs/software/geneid/Plantilla.html";
#my $FRAME2 = "$APACHE_ROOT/htdocs/software/geneid/Plantilla2.html";
my $FRAME2 = "Plantilla2.html";

my $cgi = new CGI;

#my $script = "/home/ug/gmaster/projects/astalavista/gphase3.jar";
#my $script = "~gmaster/projects/astalavista/gphase3.jar";
#my $script = "$HOME/projects/astalavista/gphase3.jar";
#my $script = "$HOME/projects/astalavista/gphase5.jar";
my $astafunk = "astafunk.jar";
my $hmmsearch= "../bin/hmmsearch";
my $hmmfetch= "../bin/hmmfetch";

#my $java = "/usr/local/Install/jdk1.5.0_01/bin/java";
#corrected for new server:
my $java = "/usr/lib/jvm/java-6-oracle/jre/bin/java";

######### VITOR WAS HERE ####################################################
my $name = ((defined($cgi->param('name'))) ? $cgi->param('name') : "" );
my $job_email = ((defined($cgi->param('job_email'))) ? $cgi->param('job_email') : "" );
my $job_id = ((defined($cgi->param('job_id'))) ? $cgi->param('job_id') : "" );

system("mkdir $job_id");



##############################################################################

#my $path_to_outdir = "/usr/local/Install/apache2/htdocs/astalavista/out";
#my $path_to_outdir = "/var/web/astalavista/out";
my $path_to_outdir = "./$job_id"; ## VITOR WAS HERE
my $outname = "landscape.html";
my $urlprefix = 'http://astalavista.sammeth.net/out/';
#my $resultfile_lifetime = 28800; # results will stay min 8H
#my $resultfile_lifetime = 86400; # 24H
my $resultfile_lifetime = 259200; # 3 days
#my $resultfile_lifetime = 30;

my $outdir = tempdir( "$path_to_outdir" );
$outdir=~/\/([^\/]+)$/ || die "internal error in outdirname!\n";
my $tempdirname = $1;

my $input="";
my @INPUT="";
# input annotation file
my $gtffile = "$outdir/astalavista.gtf";
# input ID list
my $idlistfile = "$outdir/idlist.txt";
my $idlist = "";
# url of the annotation, idlist and results files
my $gtffileweb="$urlprefix$tempdirname/astalavista.gtf";
my $outnameweb="$urlprefix$tempdirname/$outname";
my $idlistfileweb="$urlprefix$tempdirname/idlist.txt";

#my $gencodefile="/usr/local/Install/apache2/htdocs/software/astalavista/datasets/gencode/44regions_genes_CHR_coord.gtf";
#my $refseqfile="/usr/local/Install/apache2/htdocs/software/astalavista/datasets/refseq/RefSeqGenes_fromUCSC.gtf";
#my $ensemblfile="/usr/local/Install/apache2/htdocs/software/astalavista/datasets/ensembl/EnsemblGenes_fromUCSC.gtf";
#my $datasetdir="$APACHE_ROOT/htdocs/software/astalavista/datasets";
my $datasetdir="/usr/lib/cgi-bin/astalavista/datasets";
#my $gencodefile="$datasetdir/gencode/44regions_genes_CHR_coord.gtf";
#my $refseqfile="$datasetdir/refseq/RefSeqGenes_fromUCSC.gtf";
#my $ensemblfile="$datasetdir/ensembl/EnsemblGenes_fromUCSC.gtf";
#my $gencodefile="$datasetdir/human.gencode.gtf";
#my $refseqfile="$datasetdir/human.refseq.gtf";
#my $ensemblfile="$datasetdir/human.ensembl.gtf";

my $cmd="";

# get cgi parameters
my $species=     ( (defined($cgi->param('species'))    ) ? $cgi->param('species') : "" );
my $pastedinput= ( (defined($cgi->param('pastedinput'))) ? $cgi->param('pastedinput') : "" );
my $upfile=      ( (defined($cgi->param('upfile'))     ) ? $cgi->param('upfile') : "" );
my $dataset=     ( (defined($cgi->param('dataset'))    ) ? $cgi->param('dataset') : "" );
my $transcripts= ( (defined($cgi->param('transcripts'))) ? $cgi->param('transcripts') : "" );
my $events=      ( (defined($cgi->param('events'))     ) ? $cgi->param('events') : "" );
my $outputformat=( (defined($cgi->param('outputformat')))? $cgi->param('outputformat') : "" );

my $datasetfile="";

print "Content-type: text/html\n\n";

print_header();
print_title("Input");
print "AStalavista input file processing... &nbsp \n";

#####     deleting old files     #####
my $my_time = time;
opendir(DIR,"$path_to_outdir");
my @dir= sort(grep(/astalavista\.......$/,readdir(DIR)));
closedir(DIR);
foreach(@dir){
    if ((stat("$path_to_outdir/$_"))[9] <= ($my_time - $resultfile_lifetime)) {
	$cmd= "rm -R $path_to_outdir/$_";
	system($cmd);
    }
}

if($_debug){
    print STDERR "\n\n##########     DEBUGGING MODE     ##########\n\n";
}
($_debug) && print STDERR "\nFRAME2= $FRAME2\n\n";
($_debug) && print STDERR "outdir= $outdir\n";

# INPUT IS PASTED
if ($pastedinput=~ /./) {
    $input=$pastedinput;
}
# OR INPUT IS UPLOADED
elsif ($upfile =~ /./) {
    $input = $cgi->upload('upfile');
}
elsif ($dataset eq "") {
    ## ERROR: no input at all!
    print_error("<b>ERROR: Input not found..</b><br><br>Please select an input annotation (e.g. RefSeq) or paste/upload your own gtf annotation (or identifiers).");
}

open GTFFILE, ">$gtffile" or die "can't open gtf file, $gtffile!\n";
open IDLISTFILE, ">$idlistfile" or die "can't open id list file, $idlistfile!\n";

# parse the input
# special case:
goodies($input);
#@INPUT= split("[\n]",$input);
@INPUT= (($upfile =~ /./)?<$input>:split("[\n]",$input));
my $inputtype="";
my $n=0;
foreach (@INPUT) {
    $n++;
    $_ =~ s/\x0D$//; # f.ing windows/unix compatibility pb
    $_ =~ s/^[ \t\n]*//; # delete leading/trailing whitespace (spaces, tabs, empty lines)
    $_ =~ s/[ \t\n]*$//;
    if ($_=~/./) { # skip empty lines
	if ($inputtype eq "") { # first entry
	    $inputtype=(($_ =~ /^[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+/)?"gtf":"idlist");
	    ($_debug) && print "<br> line: <br> $_ <br> <br> Input type: $inputtype<br>";
	}
	if ($inputtype eq "gtf") {
	    # should be gtf format
	    ($_ =~ /^[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+/) || print_error("<b>ERROR: Problem in the GTF input:</b><br> Please check the line $n, we do not see any more the 8 tabulation-separated columns:<br>Entry= $_\n");
	    print GTFFILE "$_\n";
	}
	else { # type= idlist
	    if ($_ =~ /[ \t]/) {
		print_error("<b>ERROR: Input format error.</b><br> A gtf format requires at least 8 column (tab-separated) and a list of gene/transcript/protein identifiers requires only one entry per line<br> Please check the line $n, entry= $_\n");
	    }
	    ($_debug) && print "input: $_ <br>";
	    print IDLISTFILE "$_\n";
	}
    }
}

close GTFFILE;
close IDLISTFILE;

# at this point, we should have everything either in gtffile or in idlistfile
# in this last case, we have to extract the gtf lines from the known annotations to build the query in the gtffile.
if ($inputtype eq "idlist") {
    $cmd="cat $datasetdir/$species*gtf | grep -iwf $idlistfile > $gtffile";
    ($_debug) &&  print "command: <br> $cmd\n <br>";
    system($cmd);
    if (-z "$gtffile") {
	print_error("<b>ERROR: Gene not found.</b><br>We were unable to find any gene from the list of identifiers you provided.<br>Please make sure that you selected the right organism/genome build.");
    }
}

##### THE INPUT SHOULD CONTAIN A DATASET
my $cheatingmode=0;
if ($dataset =~ /./) { #defined ($cgi->param ('dataset')) && ($cgi->param ('dataset') =~ /./)) {
    $datasetfile="$datasetdir/"."$species"."_$dataset"."_annotation.gtf";
#    my $dataset=($cgi->param ('dataset'));
    if ((-z "$gtffile") && ($transcripts eq "all") && ($events eq "all") && ($species =~ /^human/) && ($outputformat eq "gtf")) {
	# no other input, empty gtf, and default options for human -> cheating mode
	$cheatingmode=1;
	($_debug) && print "<br> dataset on and cheating mode on<br>";
	$gtffile=$datasetfile;
	$outnameweb =~ s/$outname/\.\.\/\.\.\/datasets\/$datasetfile\.output\/$outname/;
	$outnameweb = "$urlprefix"."../datasets/$species"."_$dataset"."_annotation.gtf.output/$outname";
	$gtffileweb="$urlprefix"."../datasets/$gtffile";
	$gtffileweb =~ s/datasets.+datasets/datasets/; # marche pas?!
    }
    else {
	($_debug) && print "<br> dataset on and cheating mode off<br>";
	$cmd="cat $datasetfile >> $gtffile";
	system($cmd);
    }
}

if (-z "$gtffile") {

    if ($_debug) {
	print STDERR "Empty input file, $gtffile...\n";
    }
    print_error("<b>ERROR: GTF annotation not found..</b><br><br>Please paste or upload a GTF annotation.");
}

# print the input
#print "Your input gtf file is <a href=$urlprefix$tempdirname/astalavista.gtf>there</a><br>\n";
print "<b>completed</b><br><br>\n";
print "<b><a href=$gtffileweb>To check your input gtf annotation file, please click here.</a></b><br>\n";

my $options = "";

$options.=" -nonmd ";

if ($species =~ /./) {
#    my $species=$cgi->param ('species');
    $options.= " -genome $species";
}

if ($outputformat =~ /./) {
    $options.= " -output $outputformat ";
}

if ($transcripts eq "coding") {
    $options.=" -codingTranscripts";
}
elsif ($transcripts eq "noncoding") {
    $options.=" -noncodingTranscripts";
}
else {
    ($transcripts eq "all") || die "Inconsistent transcript option!!\n";
}

if ($events eq "coding") {
    $options.=" -CDS";
}
elsif ($events eq "noncoding") {
    $options.=" -UTR";
}
else {
    ($events eq "all") || die "Inconsistent events option!!\n";
}
$options.=" -html";
#die "\n <br> option= $options\n\n <br>";
#my $email    = $cgi->param ('email');

#my $path_to_script = "/home/ug/gmaster/projects/moby/prod/scripts/workflows_implementations";

($_debug) && print STDERR "launching the structurator...\n";

print_title("Results");
print "Alternative splicing events extraction processing... &nbsp \n";

# my $picture = qx/$path_to_script\/GenesClustering_FASTA.pl -x 1 -c $path_to_script\/workflow.config -d $matrix -t $threshold -m $method -f $gtffile >& \/dev\/null/;
#my $outmess = qx/$java -jar $script $outdir\/astalavista.gtf -html/;
my $outmess;

($_debug) && print STDERR "COMMAND LINE:\n$java -jar $script $gtffile $options\n";

eval {
#    $outmess = qx/$java -jar $script $outdir\/astalavista.gtf -html 2>&1/;
#    open(TMPOUT,">/usr/local/Install/apache2/htdocs/software/as/tmp.out")||die "can't open tmpout\n";
# /usr/local/Install/jdk1.5.0_01/bin/java  -Xmx2G -jar /home/ug/gmaster/projects/astalavista/gphase.jar
    my $cmd="$java -Xmx4G -jar $script $gtffile $options -html 2>&1";
#    print "$cmd\n";
#print "<br> <small> (command line: $cmd )</small><br>";
    ## dont compute results if already available
    if (!($cheatingmode)) {
	($_debug) && print "command: <br> $cmd\n <br>";
	$outmess = qx/$cmd/;
    }
};
#if ($@) {die $@;}
if ($@) {
    print_error("<b>ERROR: astalavista returned an error message:</b> <br><br> $@ ");
}

($_debug) && print STDERR "execution done!!\n";
if (defined $outmess && (length $outmess > 0)) {

    if ($_debug) {
	print STDERR "got a outmess!\n";
    }

##     if (!$_debug) {
##	unlink $gtffile;
##     }
##     print "Content-type: image/png\n\n";
    die $outmess;

#    print_error("<b>ERROR> astalavista retu
}
## else {
##     if (!$_debug) {
##	unlink $gtffile;
##     }
##     print "Content-type: text/html\n\n";
##     print_error("<b>ERROR> The execution of the genes clustering workflow has failed!");
## }

print "<b>completed</b><br><br>\n";
#print "<b><a href=$urlprefix$tempdirname/$outname>Here are your results</a>.<br><br>\n";
print "<b><a href=$outnameweb>To see the corresponding alternative splicing landscape, please click here</a>.<br><br>\n";
print_trailer();

#################################################

sub print_error {
# el parametro es un mensaje de error
    my ($mess) = @_;

    print "<br><center><font size=4>Sorry, an error has occured while processing your request<br><br>";
    print $mess;
    print "<br>";
    print "</FONT>";
    print "<br>";
    print 'For more information, please try our <a href="http://astalavista.sammeth.net/FAQ.html">FAQ</a> or contact the authors at: sylvain.foissac :at: crg.es';
    print "<br>";



    print '</td><td class="userspace" border=0 cellpadding=0 cellspacing=0 width="10px">&nbsp;</td></tr></table>';
#    print "</TD>";
#    print "</TR>";
#    print "</TABLE>";
#    print "<P>";
#    print "<hr>";
#    print "<p>";

    print_trailer();

    exit(1);
}

#################################################

sub goodies {
    my $input="@_";
    $input =~ s/[\n\x0D]//g;
    my $goodies=0;
    if ($input =~ /^authors$/i) {
	print "<br><center><font size=5 color=red>Congratulations!<br>You have found one of the AStalavista goodies!</font></center><br>";
	print '<br><CENTER><IMG SRC="http://astalavista.sammeth.net/pics/aut.jpg" HEIGHT=600 WIDTH=750></CENTER><br>';
	$goodies=1;
    }
    if ($input =~ /^micha$/i) {
	print "<br><center><font size=5 color=red>Congratulations!<br>You have found one of the AStalavista goodies!</font></center><br>";
	print '<br><CENTER><IMG SRC="http://astalavista.sammeth.net/pics/mic.jpg" HEIGHT=600 WIDTH=750></CENTER><br>';
	$goodies=1;
    }
    if ($input =~ /^sylvain$/i) {
	print "<br><center><font size=5 color=red>Congratulations!<br>You have found one of the AStalavista goodies!</font></center><br>";
	print '<br><CENTER><IMG SRC="http://astalavista.sammeth.net/pics/syl.jpg" HEIGHT=600 WIDTH=750></CENTER><br>';
	$goodies=1;
    }
    if ($input =~ /^roderic$/i) {
	print "<br><center><font size=5 color=red>Congratulations!<br>You have found one of the AStalavista goodies!</font></center><br>";
	print '<br><CENTER><IMG SRC="http://astalavista.sammeth.net/pics/rod.jpg" HEIGHT=600 WIDTH=750></CENTER><br>';
	$goodies=1;
    }
    if ($input =~ /^miss[\s-]splicing$/i) {
	print "<br><center><font size=5 color=red>Congratulations!<br>You have found one of the AStalavista goodies!</font></center><br>";
	print '<br><CENTER><IMG SRC="http://astalavista.sammeth.net/pics/mis.jpg" HEIGHT=600 WIDTH=750></CENTER><br>';
	$goodies=1;
    }
    if ($input =~ /^smiley$/i) {
	print "<br><center><font size=5 color=red>Congratulations!<br>You have found one of the AStalavista goodies!</font></center><br>";
	print '<br><CENTER><IMG SRC="http://astalavista.sammeth.net/pics/smi.gif" HEIGHT=600 WIDTH=750></CENTER><br>';
	$goodies=1;
    }
    if ($goodies==1) {
	print "</FONT>";
	print "</TD>";
	print "</TR>";
	print "</TABLE>";
	print "<P>";
	print "<hr>";
	print "<p>";
	print_trailer();
	exit(1);
    }
}

#################################################

sub print_header {
my $HEADER =<<'+++EOH+++';
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
<html>
<head>
 <META HTTP-EQUIV="Content-Type" CONTENT="text/html; charset=iso-8859-1">
 <META NAME="Generator" CONTENT="Hacker Hands -from geneid server">
 <META NAME="author" CONTENT="Michael Sammeth and Sylvain Foissac">
 <!-- KEYWORDS and DESCRIPTION for INTERNET SEARCH ENGINES -->
  <META NAME="keywords"    CONTENT="alternative splicing, astalavista server, michael sammeth, roderic guigo, sylvain foissac">
  <META NAME="description" CONTENT="AStalavista web server.">
 <!-- KEYWORDS and DESCRIPTION for INTERNET SEARCH ENGINES -->
 <TITLE>AStalavista - Alternative Splicing landscape</TITLE>
<!-- <BASE HREF="/"> -->
<link rel=stylesheet type="text/css" href="/Genome.css" title="Genome">
<STYLE type="text/css">
 @import url(/Genome_forms_mozilla.css);
</STYLE>
<SCRIPT type="text/javascript" src="/Genome.js"></SCRIPT>
</head>

<body onResize="init();" MARGINWIDTH="2" MARGINHEIGHT="0">
<A NAME="TOC"></A>

<!-- ############### TITLE AREA ############### -->
<table border=0 cellpadding=0 cellspacing=0 width="100%">
<tr>
<td align="center">
<FONT size=4 class="tit">Genome BioInformatics Research Lab</FONT>
</td>
</tr>
</table>

<!-- ############### MENU AREA ############### -->
<table class="default" border=0 cellpadding=0 cellspacing=0 width="100%">
<tr>
<td class="links" align="CENTER">
 <a href="http://genome.imim.es/main/help.html"
  onMouseover="window.status='Asking for Help';">Help</a>
&nbsp;<FONT style="color: #005050;">&#124;</FONT>&nbsp;
 <a href="http://genome.imim.es/main/news.html"
  onMouseover="window.status='News';">News</a>
&nbsp;<FONT style="color: #005050;">&#124;</FONT>&nbsp;
 <a href="http://genome.imim.es/main/people.html"
  onMouseover="window.status='Who is who in Genome Informatics Research Laboratory';">People</a>
&nbsp;<FONT style="color: #005050;">&#124;</FONT>&nbsp;
 <a href="http://genome.imim.es/main/research.html"
  onMouseover="window.status='Research at our Group';">Research</a>
&nbsp;<FONT style="background: #005050; color: #FFFFFF;">&nbsp;
 <a href="http://genome.imim.es/main/software.html"
  style="background: #005050; color: #DDDDDD;"
  onMouseover="window.status='Software developed in our Group';">Software</a>
&nbsp;</FONT>&nbsp; <!-- Current menu -->
 <a href="http://genome.imim.es/main/publications.html"
  onMouseover="window.status='Publications by our Group';">Publications</a>
&nbsp;<FONT style="color: #005050;">&#124;</FONT>&nbsp;
 <a href="http://genome.imim.es/main/links.html"
  onMouseover="window.status='Interesting Links';">Links</a>
</td>
</tr>
<tr>
<td class="links" align="CENTER">
 <a href="http://genome.imim.es/main/databases.html"
  onMouseover="window.status='Datasets created/curated by our Group';" class="hh">Resources &amp; Datasets</a>
&nbsp;<FONT style="color: #005050;">&#124;</FONT>&nbsp;
 <a href="http://genome.imim.es/genepredictions/index.html"
  onMouseover="window.status='Gene Predictions produced by our Software';" class="hh">Gene Predictions</a>
&nbsp;<FONT style="color: #005050;">&#124;</FONT>&nbsp;
 <a href="http://genome.imim.es/main/seminars.html"
  onMouseover="window.status='Seminars, Courses and Group Sessions';" class="hh">Seminars &amp; Courses</a>
</td>
</tr>
</table><BR class="hhh">

<!-- ############### PATH AREA ############### -->
<table border=1 cellpadding=0 cellspacing=0 width="100%">
<tr>
<td class="path" border=0 cellpadding=0 cellspacing=0 ALIGN="LEFT">
 <table border=0 cellpadding=0 cellspacing=1>
 <tr>
 <td class="ppath" border=0 cellpadding=0 cellspacing=0>
&nbsp;
<a href="http://www.imim.es/"
 onMouseover="window.status='Institut Municipal d\'Investigaci&oacute; M&egrave;dica';">IMIM</a>
 </td><td class="ppath" border=0 cellpadding=0 cellspacing=0>
<IMG SRC="http://genome.imim.es/g_icons/dmd.gif" ALT="*" HEIGHT=8 WIDTH=8 BORDER=0>
 </td><td class="ppath" border=0 cellpadding=0 cellspacing=0>
<a href="http://www.upf.edu/"
 onMouseover="window.status='Universitat Pompeu Fabra';">UPF</a>
 </td><td class="ppath" border=0 cellpadding=0 cellspacing=0>
<IMG SRC="http://genome.imim.es/g_icons/dmd.gif" ALT="*" HEIGHT=8 WIDTH=8 BORDER=0>
 </td><td class="ppath" border=0 cellpadding=0 cellspacing=0>
<a href="http://www.crg.es/"
 onMouseover="window.status='Centre de Regulaci&oacute; Gen&ograve;mica';">CRG</a>
 </td><td class="ppath" border=0 cellpadding=0 cellspacing=0>
<IMG SRC="http://genome.imim.es/g_icons/dmd.gif" ALT="*" HEIGHT=8 WIDTH=8 BORDER=0>
 </td><td class="ppath" border=0 cellpadding=0 cellspacing=0>
<a href="http://www.imim.es/grib/eng/default.htm"
 onMouseover="window.status='Grup de Recerca en InformÃ£Â tica Biom&egrave;dica';">GRIB</a>
 </td><td class="ppath" border=0 cellpadding=0 cellspacing=0>
<A HREF="http://genome.imim.es/index.html" onMouseover="window.status='Genome Informatics Research Laboratory HOME PAGE: You are Welcome!!!'">
<IMG SRC="http://genome.imim.es/g_icons/lhr.gif" ALT="HOME" HEIGHT=10 WIDTH=36 BORDER=0>
</A>
 </td><td class="ppath" border=0 cellpadding=0 cellspacing=0>
 <a href="http://genome.imim.es/main/software.html"
  onMouseover="window.status='Software developed in our Group';">Software</a>
 </td>
 <td class="ppath" border=0 cellpadding=0 cellspacing=0>
<IMG SRC="http://genome.imim.es/g_icons/dmd.gif" HEIGHT=8 WIDTH=8 BORDER=0>
</td><td class="ppath" border=0 cellpadding=0 cellspacing=0>
<a href="http://astalavista.sammeth.net" onMouseover="window.status='AStalavista server';">AStalavista</a>
</td>

 </tr>
 </table>
</td>
</tr>
</table>

<!-- ############### USERS AREA ############### -->
<table border=0 cellpadding=0 cellspacing=0 width="100%">
<tr border=0 cellpadding=0 cellspacing=0>
<td class="section" border=0 cellpadding=0 cellspacing=0 width="10px">&nbsp;</td>

<td class="section" ALIGN=CENTER>
<CENTER><FONT size=6 class="tgen">- AStalavista web server -<br>Alternative Splicing transcriptional landscape visualization tool</FONT></CENTER>
</td>
<td class="section" border=0 cellpadding=0 cellspacing=0 width="10px">&nbsp;</td>
</tr>
<tr border=0 cellpadding=0 cellspacing=0>
<td class="userspace" border=0 cellpadding=0 cellspacing=0 width="10px">&nbsp;</td>
<td class="userspace" border=0 cellpadding=0 cellspacing=0>

<!-- <FORM method="POST" enctype="multipart/form-data" action="/cgi-bin/as/astalavista.cgi">-->

<!--<br><CENTER><FONT size=2><b>From a genomic annotation file, identify and display the set of alternative splicing events.</b></FONT></CENTER><br> -->
+++EOH+++

print "$HEADER";
}

sub print_title {
my ($title) = @_;
my $TITLE1 =<<'+++EOH+++';
<br>
<TABLE border=0 cellpadding=0 cellspacing=0 width=100%>
<TR>
<TD class="section" border=0 cellpadding=0 cellspacing=0 width=80% ALIGN=LEFT>
<FONT size=5 class="K">
+++EOH+++

my $TITLE2 =<<'+++EOH+++';
</FONT>
</TD>

<TD class="section" border=0 cellpadding=0 cellspacing=0 width=20% ALIGN=RIGHT>
<a href="index.html#TOP"
 onMouseover="window.status='TOP:INDEX';">
<IMG class="pnt" SRC="/g_icons/top.gif" HEIGHT=15 WIDTH=15 BORDER=0></a>
</TD></TR>
</TABLE>
<br>
+++EOH+++

print "$TITLE1";
print "$title";
print "$TITLE2";
}

sub print_trailer {
my $HEADER =<<'+++EOH+++';
<br>
</td>
<td class="userspace" border=0 cellpadding=0 cellspacing=0 width="10px">&nbsp;</td>
</tr>
</table>
<!-- ############### TRAILER AREA ############### -->
<table class="trailer" border=0 cellpadding=0 cellspacing=0 width="100%">
 <tr class="default">
 <td id="traila" align=left>
&nbsp;
<A HREF="http://genome.imim.es/main/disclaimer.html"
 onMouseover="window.status='Data Release Policy and Software License';">Disclaimer</A>
 </td>
 <td id="trailb" align=center>
<SCRIPT type="text/javascript" languaje="javascript">
<!--
 LastUpdateSplit()
// -->
</SCRIPT>
 </td>
 <td id="trailc" align=rigth>
<A HREF="MAILTO:gmaster@imim.es"
 onMouseover="window.status='Contact with our webmaster';">webmaster</A>
&nbsp;
 </td>
 </tr>
</table>

<SCRIPT type="text/javascript" languaje="javascript">
<!--
 reportErrors();
// -->
</SCRIPT>
<!-- ############### EOF ############### -->
</body>
</html>
+++EOH+++
print "$HEADER";
}
