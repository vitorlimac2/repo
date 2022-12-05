#!/usr/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
#use CGI qw/:push -nph/;
#use CGI qw(:standard -nph);
use Data::Dumper; # for debugging input to cgi
use File::Copy;
#use File::Temp qw/tempfile/;
use File::Temp qw/tempdir/;

my $_debug = 1;
#warn Dumper [<STDIN>]; # clears CGI variables

#my $APACHE_ROOT = $ENV{'APACHE_ROOT'};
#my $HOME = "~gmaster";
# $ENV{'HOME'};

my $script = "../tools/astafunk.gui/barna.gui.jar";
my $astafunk = "../tools/astafunk.jar";
my $hmmsearch= "../tools/hmmsearch";

my $path_to_outdir = "../jobs";
my $outname = "index.html";
my $urlprefix = 'http://afunk.sammeth.net/jobs/';
#my $resultfile_lifetime = 259200; # 3 days
my $resultfile_lifetime = 180;

## setup temp dir for job
my $outdir = tempdir( "$path_to_outdir/astafunk.XXXXXX" );
$outdir=~/\/([^\/]+)$/ || die "internal error in outdirname!\n";
my $tempdirname = $1;
chmod 0755, $outdir;
my $ppage = "progress.html";
my $pfile = "$outdir/$ppage";
my $purl = "$urlprefix$tempdirname/$ppage";

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

my $datasetdir="../data";
my $testoutfile="afunk_all_output.txt";

my $cmd="";

# get cgi parameters
#my $species=     ( (defined($cgi->param('species'))    ) ? $cgi->param('species') : "" );

#my $datasetfile="../data/afunk_all_output.txt";

my $options="";
my $outmess="";
my $fh;
my $cheatingmode=0;


## DO CGI STUFF

print "Content-type: text/html\n\n";

my $cgi = new CGI;
#CGI->nph(1); # NPH outdated.. use dynamically updated progress html

#my $hmm_file = $cgi->param("hmmToUpload");
#$hmm_file =~m/^.*(\\|\/)(.*)/;
#write_file($hmm_file);
#my $uf = $cgi->param('upfile');
my $upname = $cgi->param('file_upload');
my $upfile = $cgi->upload('file_upload');
my $some_file = "$outdir/$upname";
print "upload $upname , $upfile to $some_file\n";
upload($upfile,$some_file);


# GENOME
my $genome_file = $cgi->param("organismToUpload");
print "$genome_file";
#$genome_file =~m/^.*(\\|\/)(.*)/;
#write_file($genome_file);

# TRANSCRIPTOME
#my $gtf_file = $cgi->param("annotationToUpload");
#$gtf_file =~m/^.*(\\|\/)(.*)/;
#write_file($gtf_file);

# HMM
#my $hmm_file = $cgi->param("hmmToUpload");
#$hmm_file =~m/^.*(\\|\/)(.*)/;
#write_file($hmm_file);


## PROGRESS PAGE
open($fh, ">$pfile") or die "Cannot open $pfile.";
print $fh "<HTML>\n<HEAD>\n<META http-equiv=\"refresh\" content=\"3\">\n";
print $fh "<META http-equiv=\"Pragma\" content=\"no-cache\">\n</HEAD>";
print $fh "<BODY>\n<h1>Astafunk in Progress..</h1>\n";
print $fh "<P>Depending on the input you provided, this is likely to take a while. ";
print $fh "A link to the results will appear at the end of page once they are ready.</P>\n";
close $fh or warn "Something went wrong with writing to $pfile.";

print "Content-type: text/html\n\n";
print "<H1>Job Submitted</H1>\n";
print "<BR><BR>";
print "<I>Your job has been sucessfully submitted and can be monitored via the link:</I><BR><BR>\n";
print "&nbsp;&nbsp;&nbsp;<A href=\"$purl\">$purl</A><BR><BR>\n";
print "<I>At your consideration bookmark the page above in case of long running processes.<BR><BR>";

my $genome_dir;
my $ref_tx_fasta = "$outdir/reference_tx.fasta";
my $ref_domain_tab = "$outdir/reference_file";
my $afunk_output = "$outdir/output.txt";

## FORK OFF CHILD
if (fork() == 0) {

    # From here on, we're in the child..
    # close stdout/stderr so browser doesn't wait on child
    close STDOUT;
    close STDERR;

    # run stuff and append to progress html
    my $fh;
    open($fh, ">>$pfile"); # or..
    print $fh "<H1>Running Astafunk..</H1>\n";
    print $fh "<I>Detecting AS events and protein domains in your input. ";
    print $fh "This is likely to take long.</I><BR>\n";
    close $fh; # or..

    # PIPELINE

    # print reference transcripts
    #my $cmd = "java -Xmx4G -jar $astafunk --astref --gtf $gtf_file --genome $genome_dir > $ref_tx_fasta";
    #poll($cmd,$pfile);

#my $cmd1_to_print_reference_transcripts = "$java -Xmx4G -jar $astafunk --astref --gtf $gtf --genome $genome_dir > $ref_tx_fasta";
#my $cmd2_to_create_ref_table = "$hmmsearch --cpu 10 --noali --cut_ga --domtblout $ref_domain_tab $pfam $ref_tx_fasta";
#my $cmd3_to_run_afunk = "$java -Xmx4G -jar $astafunk --cpu 10 --all --reference --gtf $gtf --genome $genome_dir > $afunk_output";

    copy("$datasetdir/$testoutfile",$outdir) or print $fh "<B>Copy failed: $!</B><BR>";
    sleep(1); # :)
    my $resultweb = "$urlprefix$tempdirname/$testoutfile";

    open($fh, ">>$pfile"); # or..
    print $fh "<H2>..done!</H2>";
    print $fh "<H1>Formatting Results..</H1>\n";
    close $fh; # or..

    # start process and poll output
    my $cmd="java -jar $script $outdir/$testoutfile 2>&1";
    poll($cmd,$pfile);

    open($fh, ">>$pfile"); # or..
    #print $fh $outmess;
    print $fh "<H2>..done!</H2>\n";
    print $fh "<H1>Ready.</H1>\n";
    print $fh "<I>Your results are availalble at the page:</I><BR><BR>\n";
    print $fh "&nbsp;&nbsp;&nbsp;<A href=$outnameweb>$outnameweb</A><BR><BR>\n";
    print $fh "<I>or download a tab-delimited file with the results:</I><BR><BR>\n";
    print $fh "&nbsp;&nbsp;&nbsp;<A href=$resultweb>$resultweb</A><BR><BR>\n";
    print $fh "<I>We keep your files for more or less $resultfile_lifetime seconds on our servers, ";
    print $fh "please make sure to take out important results before that time has passed!</I><BR>\n";
    print $fh "<I>ASta!</A>\n";
    print $fh "</BODY>\n</HTML>";

    close $fh; # or..
    exit;  # ends the child process.
}

## CGI TERMINATES HERE
sleep(3); # wait 3 seconds to allow the child to start up
print "Asta for now..</I><BR>";
exit; # parent (cgi) exits here


### SUBROUTINES

sub write_file{
    my ($file) = @_;
    open ( LOCAL, ">$outdir/$file" ) or die "$!";
    binmode LOCAL;

    while ( <$file> )
    {
        print LOCAL $_;
    }

    close LOCAL;

    print "$file has been successfully uploaded... thank you.\n";
}

sub poll{
    my $cmd = shift;
    my $pfile = shift;
    open my $cmd_fh, "$cmd |";
    while (<$cmd_fh>) {
    	open($fh, ">>$pfile"); # or..
    	print $fh "$_\n";
    	close $fh; # or..
     	sleep(1);
    }
    close $cmd_fh;
}

sub upload{
    my $io_handle = shift;
    my $some_file_path = shift;
    print "copying to $some_file_path\n";
    #my $buffer;
    open ( my $out_file, ">$some_file_path" ) or die $!;
    #while ( my $bytesread = $io_handle->read(my $buffer,1024) ) {
    #	print $out_file $buffer;
    #}
    while (read ($io_handle, my $Buffer, 1024)) {
    	print $out_file $Buffer;
    }
    close $out_file;

}

sub get_file{
    my $File = shift; #$req->param('file');
    my $outdir = shift;
    $File =~ /.*\."?(\w*)"?$/;
    my $Filename = lc($1);

    open (FILE, ">$outdir/$Filename");
    while (read ($File, my $Buffer, 1024)) {
        print FILE $Buffer;
    }
    close FILE;
}
  
