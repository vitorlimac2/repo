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
#use IO::Zlib;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
#use File::MMagic;
#use Archive::Tar;
use Compress::Zlib;

my $_debug = 1;
#warn Dumper [<STDIN>]; # clears CGI variables
my $cmd="";
my $outmess="";
my $fh;
my $cheatingmode=0;

#my $APACHE_ROOT = $ENV{'APACHE_ROOT'};
#my $HOME = "~gmaster";
# $ENV{'HOME'};

# cgi field names
my $genome_sel_id = "genome_sel";
my $genome_ul_id = "genome_ul";
my $trptome_sel_id = "trptome_sel";
my $trptome_ul_id = "trptome_ul";
my $domains_sel_id = "domains_sel";
my $domains_ul_id = "domains_ul";
my $gene_ids = "gene_ids";
my $domain_ids = "domain_ids";

# server directories and files
my $jobsdir = "../jobs";
my $datadir = "../data";

my $f_genome = "genome.mfasta";
my $f_trptome = "trptome.gtf";
my $f_domains = "domains.hmm";
my $f_domdesc = "domdesc.tsv";
my $f_query_genes = "query_genes.txt";
my $f_query_domains = "query_domains.txt";
my $f_gids = "ids_genes.txt";
my $f_dids = "ids_domains.txt";
my $testoutfile="afunk_all_output.txt";

my $outname = "index.html";
my $urlprefix = 'http://afunk.sammeth.net/jobs/';
#my $resultfile_lifetime = 259200; # 3 days
my $resultfile_lifetime = 180;

my $script = "../tools/astafunk.gui/barna.gui.jar";
my $astafunk = "../tools/astafunk.jar";
my $hmmsearch= "../tools/hmmsearch";
my $gtf_printGidTid = "./gtf_print_gid_tid.sh";
my $filter_hmm = "./filterHMM.sh";

## DO CGI STUFF

print "Content-type: text/html\n\n";

my $time_t0= time();
my $timestr= localtime($time_t0);
print "<H2>Gathering Input Data</H2>\n";
print "<I>Request received:</I>&nbsp;<font size=\"-1\">$timestr</font><BR>\n";
my $outdir = tempdir( "$jobsdir/astafunk.XXXXXX" );
#my $outdir_job=~/\/([^\/]+)$/ || die "<B>Error setting up output folder!<\B><BR>\n";
my @outdir_tokens = split /\//, $outdir || die "<B>Error setting up output folder!<\B><BR>\n";
my $outdir_name = $outdir_tokens[2];
#chmod 0755, $outdir; # DEBUG only
my $outdir_url = "$urlprefix$outdir_name";
my $p_genome = "$outdir/$f_genome";
my $p_trptome = "$outdir/$f_trptome";
my $p_domains = "$outdir/$f_domains";
my $p_domdesc = "$outdir/$f_domdesc";
my $p_tmp = "$outdir/temp.tmp";
my $p_qgene = "$outdir/$f_query_genes";
my $p_qdom = "$outdir/$f_query_domains";
my $p_gids = "$outdir/$f_gids";
my $p_dids = "$outdir/$f_dids";

print "<UL>";
print "<LI><I>Set up the working directory for job:</I>&nbsp;$outdir_name</LI>";
my $cgi = new CGI; # CGI is uploading right here already
#CGI->nph(1); # NPH outdated.. use dynamically updated progress html
my $time_t1 = ptime("<LI><I>Uploaded data.</I>", $time_t0, \*STDOUT);
my $genome_sel = $cgi->param($genome_sel_id);
my $genome_ul = defined($cgi->upload($genome_ul_id)) ? $cgi->upload($genome_ul_id) : "";
my $b_genome = ($genome_ul eq "") ? cdata("$datadir/$genome_sel",$p_genome) : upload($genome_ul,$p_genome);
print (($genome_ul eq "") ? "<LI><I>Genome ready:</I>&nbsp$genome_sel</LI>" : "<LI><I>Genome uploaded:</I>&nbsp;$genome_sel</LI>");
my $trptome_sel = $cgi->param($trptome_sel_id);
my $trptome_ul = defined($cgi->upload($trptome_ul_id)) ? $cgi->upload($trptome_ul_id) : "";
my $b_trptome = ($trptome_ul eq "") ? cdata("$datadir/$trptome_sel",$p_trptome) : upload($trptome_ul,$p_trptome);
print (($trptome_ul eq "") ? "<LI><I>Transcriptome ready:</I>&nbsp;$trptome_sel</LI>" : "<LI><I>Transcriptome uploaded:</I>&nbsp;$trptome_ul</LI>");
my $domains_sel = $cgi->param($domains_sel_id);
my $domains_ul = defined($cgi->upload($domains_ul_id)) ? $cgi->upload($domains_ul_id) : "";
my $b_domains = ($domains_ul eq "") ? cdata("$datadir/$domains_sel",$p_domains) : upload($domains_ul,$p_domains);
print (($domains_ul eq "") ? "<LI><I>Domains ready:</I>&nbsp;$domains_sel</LI>" : "<LI><I>Domains uploaded:</I>&nbsp;$domains_ul</LI>");

my $gids= ( (defined($cgi->param($gene_ids))) ? $cgi->param($gene_ids) : "" );
if ($gids ne "") { 
    my $n = printIDs($gids, $p_qgene, \*STDOUT);
    print "<LI><I>Gene filters:</I>&nbsp;received $n queries</LI>";
}
my $dids= ( (defined($cgi->param($domain_ids))) ? $cgi->param($domain_ids) : "" );
if ($dids ne "") { 
    my $n = printIDs($dids, $p_qdom, \*STDOUT);
    print "<LI><I>Domain filters:</I>&nbsp;received $n queries</LI>";
}

print "</UL><BR>\n";
if ($_debug) {
   print "Start: time=$time_t0<BR>\n";
   print "Genome sel=$genome_sel, ul=$genome_ul, dest=$p_genome<BR>\n";
   print "Transcriptome sel=$trptome_sel, ul=$trptome_ul, dest=$p_trptome<BR>\n";
   print "Domains sel=$domains_sel, ul=$domains_ul, dest=$p_domains<BR>\n";
}


#exit;

## CHECK STUFF


## PROGRESS PAGE

my $ppage = "progress.html";
my $pfile = "$outdir/$ppage";
my $purl = "$outdir_url/$ppage";

open($fh, ">$pfile") or die "Cannot open $pfile.";
print $fh "<HTML>\n<HEAD>\n<META http-equiv=\"refresh\" content=\"3\">\n";
print $fh "<META http-equiv=\"Pragma\" content=\"no-cache\">\n</HEAD>";

print $fh "<BODY>\n<h1>Job $outdir_name</h1>\n";
print $fh "<P><I>Depending on the input you provided, AStafunk may take a while to complete. This page will update dynamically and ";
print $fh "eventually links to the results will appear at the end of the page once all computations terminated successfully.</I></P>\n";
close $fh or warn "Something went wrong with writing to $pfile.";

#print "Content-type: text/html\n\n";
print "<H2>Job Submitted!</H2>\n";
#print "<BR><BR>";
print "<I>Your job</I> &quot;$outdir_name&quot; <I>has been sucessfully submitted and can be monitored on the page:</I><BR><BR>\n";
print "&nbsp;&nbsp;&nbsp;<A href=\"$purl\">$purl</A><BR><BR>\n";
print "<I>At your convenience, bookmark the link above in order to track the potentially long running process.<BR><BR>";

my $genome_dir = $outdir;
my $ref_tx_fasta = "$outdir/reference_tx.fasta";
my $ref_domain_tab = "$outdir/reference_file";
my $afunk_output = "$outdir/output.txt";

## FORK OFF CHILD
if (fork() == 0) {

    # From here on, we're in the child..
    # close stdout/stderr so browser doesn't wait on child
    close STDOUT;
    close STDERR;
    my $fh;

    # Preprocess files
    open($fh, ">>$pfile"); # or..
    print $fh "<H2>Preprocessing</H2>\n";
    print $fh "<I>Preparing input files for AStafunk run.</I><BR>\n";
    print $fh "<UL>\n";

#    print $fh "<LI><I>Checking genome:</I>&nbsp;";
#    my $cmd = "file $p_genome";
#    my $res = qx/$cmd/;
#    if ($_debug) {
#	print $fh "<BR>$cmd: $res<BR>";
#    }
#    close $fh; # or..
#    open($fh, ">>$pfile"); # or..
#    if ($res =~ m/gzip compressed data/) {
#	$time_t1 = time();
#	gunzip $p_genome => $p_tmp or return undef;
#	rename $p_tmp => $p_genome;
#	ptime("gunzipped. ", $time_t1, $fh);
#    } else {
#	print $fh "not gzip compressed.";
#    }
#    print $fh "</LI>\n";
#    close $fh; # or..

    # genome
    inflate("Inflating genome", $p_genome, $p_tmp, $pfile);
    open($fh, ">>$pfile"); # or..
    print $fh "<LI><I>Splitting genome:</I>&nbsp;";
    $time_t1 = time();
    my $cmd = "../tools/separate_multifasta.sh $outdir/$f_genome $outdir";
    my $res = qx/$cmd/;
    ptime("found $res contigs. ", $time_t1, $fh);
    print $fh "</LI>\n";
    close $fh; 

    # trptome
    inflate("Inflating transcriptome", $p_trptome, $p_tmp, $pfile);
    if ($gids ne "") {
    	$cmd = "grep --ignore-case -f $p_qgene $p_trptome | $gtf_printGidTid 2>&1 | sort | uniq >$p_gids; wc -l $p_gids | awk '{print \$1}'";
    	$res = qx/$cmd/;
    	pprint($pfile, "<LI><I>Filtering transcriptome:</I>&nbsp;queries matched $res gene and transcript IDs, ");
    	if ($_debug) {
   	    pprint($pfile, "<BR><B>Command:&nbsp;</B>$cmd<BR>");
    	}
    	$cmd = "grep -f $p_gids $p_trptome | awk '\$3==\"exon\" || \$3==\"CDS\"' > $p_tmp; wc -l $p_tmp | awk '{print \$1}'";
    	if ($_debug) {
            pprint($pfile, "<BR><B>Command:&nbsp;</B>$cmd<BR>");
    	}
    	$res = qx/$cmd/;
    	pprint($pfile, "$res GTF lines.</LI>");
    	rename $p_tmp => $p_trptome;
    }

    inflate("Inflating domains", $p_domains, $p_tmp, $pfile);
    if ($dids ne "") {
        if ($domains_ul eq "") {
	    (my $basename = $domains_sel) =~ s/\.hmm.*$//;   # /\.[^.]+$//;
	    cdata("$datadir/${basename}_info.tsv",$p_domdesc);
	    if ($_debug) {
		pprint($pfile, "<BR><B>Copied:&nbsp;</B>$datadir/${basename}_info.tsv&nbsp;&gt;&nbsp;$p_domdesc<BR>");
	    }
	} else {
	    $cmd = "";
	}
        $cmd = "grep --ignore-case -f $p_qdom $p_domdesc | awk '{print \$1}' >$p_dids; wc -l $p_dids | awk '{print \$1}'";
        $res = qx/$cmd/;
        pprint($pfile, "<LI><I>Filtering domains:</I>&nbsp;queries matched $res domain profiles, ");
        if ($_debug) {
            pprint($pfile, "<BR><B>Command:&nbsp;</B>$cmd<BR>");
        }
        $cmd = "cat $p_domains | $filter_hmm $p_dids > $p_tmp; wc -l $p_tmp | awk '{print \$1}'";
        if ($_debug) {
            pprint($pfile, "<BR><B>Command:&nbsp;</B>$cmd<BR>");
        }
        $res = qx/$cmd/;
        pprint($pfile, "$res HMM lines.</LI>");
        rename $p_tmp => $p_domains;
    }

    # PIPELINE

    open($fh, ">>$pfile"); # or..
    print $fh "</UL>\n";    
    print $fh "<H2>Running Astafunk</H2>\n";
    print $fh "<UL>\n";
    close $fh;

    # reference transcripts
    my $f_reftx = "ref_transcripts.fasta";
    $cmd = "java -Xmx4G -jar $astafunk --astref --gtf $outdir/$f_trptome --genome $outdir > $outdir/$f_reftx";
    open($fh, ">>$pfile"); # or..
    print $fh "<LI><I>Extracting reference transcript sequences.</I><BR>";
    if ($_debug) {
    	print $fh "Command: $cmd</LI><BR>\n";
    }
    close $fh; # or..
    $time_t1 = time();
    poll($cmd,$pfile);
    open($fh, ">>$pfile"); # or..
    ptime("<BR>done. ", $time_t1, $fh);
    print $fh "</LI>\n";
    close $fh;

    # referene domains
    my $f_refdom = "ref_domains.txt";
    $cmd = "$hmmsearch --cpu 10 --noali --cut_ga --domtblout $outdir/$f_refdom $outdir/$f_domains $outdir/$f_reftx";
    open($fh, ">>$pfile"); # or..
    print $fh "<BR><LI><I>Collecting reference domains.</I><BR>\n";
    if ($_debug) {
    	print $fh "Command: $cmd</LI><BR>\n";
    }
    close $fh; # or..
    $time_t1 = time();
    #devnull($cmd,$pfile); # hmmer doesnt like polling
    if ($_debug) {
    	$outmess = qx/$cmd/;
    }
    open($fh, ">>$pfile"); # or..
    ptime("<BR>done. ", $time_t1, $fh);
    print $fh "</LI>\n";
    close $fh;

    # astafunk
    my $f_afunk_out = "astafunk_results.txt";
    open(my $fix, '>', "$outdir/$f_afunk_out");    
    print $fix "#";
    close $fix;
    $cmd = "java -Xmx4G -jar $astafunk --cpu 10 --all --reference $outdir/$f_refdom --hmm $outdir/$f_domains --gtf $outdir/$f_trptome --genome $outdir >> $outdir/$f_afunk_out";
    open($fh, ">>$pfile"); # or..
    print $fh "<BR><LI><I>Searching domains in AS events. Likely to take long.</I><BR>\n";
    if ($_debug) {
    	print $fh "Command: $cmd</LI><BR>\n";
    }
    close $fh; # or..
    $time_t1 = time();
    poll($cmd,$pfile);
    open($fh, ">>$pfile"); # or..
    ptime("<BR>done. ", $time_t1, $fh);
    print $fh "</LI>\n";
    close $fh;

    #copy("$datadir/$testoutfile",$outdir) or print $fh "<B>Copy failed: $!</B><BR>";
    #sleep(1); # :)

    # spliceosigner
    my $cmd="java -jar $script $outdir/$f_afunk_out 2>&1";
    open($fh, ">>$pfile"); # or..
    print $fh "<BR><LI><I>Rendering HTML visualization from results.</I><BR>\n";
    close $fh; # or..
    $time_t1 = time();
    poll($cmd,$pfile);
    open($fh, ">>$pfile"); # or..
    ptime("<BR>done. ", $time_t1, $fh);
    print $fh "</LI>\n";
    close $fh;

    # results
    my $outnameweb = "$outdir_url/index.html";
    my $resultweb = "$urlprefix$outdir_name/astafunk_results.txt";
    open($fh, ">>$pfile"); # or..
    #print $fh $outmess;
    print $fh "</UL><H2>Results</H2>\n";
    print $fh "<I>Your results are availalble on the html page:</I><BR><BR>\n";
    print $fh "&nbsp;&nbsp;&nbsp;<A href=$outnameweb>$outnameweb</A><BR><BR>\n";
    print $fh "<I>or download a tab-delimited file with the results:</I><BR><BR>\n";
    print $fh "&nbsp;&nbsp;&nbsp;<A href=$resultweb>$resultweb</A><BR><BR>\n";
    print $fh "<I>We keep your files for more or less $resultfile_lifetime seconds on our servers, ";
    print $fh "please make sure to take out important results before that time has passed!</I><BR><BR>\n";
    print $fh "<I>AStalavista!</A>\n";
    print $fh "</BODY>\n</HTML>";
    close $fh; # or..

    exit;  # ends the child process.
}

## PARENT CGI TERMINATES HERE
sleep(3); # wait 3 seconds to allow the child to start up
$timestr = localtime();
print "Exiting now:</I>&nbsp;<font size=\"-1\">$timestr</font><BR><BR>\n";
print "Astalavista!<BR>\n";
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

sub devnull{
    my $cmd = shift;
    my $pfile = shift;
    open my $cmd_fh, "$cmd |";
    my $counter = 0;
    while (<$cmd_fh>) {
        open($fh, ">>$pfile"); # or..
        my $s = "$_";
        if (length($s) > 0) {
            #if($counter % 100 == 0) {
	    #	print $fh "+";
	    #} elif($counter % 10 == 0) {
            #	print $fh "|";
	    #} else {
            #    print $fh ".";
	    #}
            $counter = $counter + 1;
            if ($counter % 300 == 0) {
		print $fh "<BR>\n";
	    }
	}
        close $fh; # or..
        sleep(5);
    }
    close $cmd_fh;
}

sub upload{
    my $src = shift;  # handle
    my $dest = shift; # path
    open ( my $f_dest, ">$dest" ) or die "upload failed! $!";
    while (read ($src, my $buffer, 1024)) {
        print $f_dest $buffer;
    }
    close $f_dest;
    return $dest;
}

sub cdata {
   my $src = shift;  # path
   my $dest = shift; # path
   copy($src,$dest) or die("<B>Copy $src > $dest failed!</B><BR> $!");
   return $dest;
}

sub tdiff {
    my( $seconds ) =  (@_ );

    if ( $seconds < 60 ) {
        # less than a minute
        return "$seconds sec.";
    }
    if ( $seconds <= ( 60 * 60 ) ) {
        # less than an hour
        my $v= int($seconds/ 60 );
        my $r= tdiff($seconds % 60);
        return( "$v min, $r" );
    }
    if ( $seconds <= ( 60 * 60 * 24 ) ) {
        # less than a day
        my $v = int( $seconds/(60 * 60) );
        my $r = tdiff($seconds % (60 * 60) );
	return( "$v hours, $r" );
    }
    if ( $seconds <= ( 60 * 60 * 24 * 7 ) ) {
        # less than a week
        my $v = int( $seconds/(60*60*24));
	my $r = tdiff($seconds % (60*60*24));
        return( "$v days, $r" );
    }

    # fall-back weeks ago
    my $v = int( $seconds/(60*60*24*7));
    my $r = tdiff($seconds % (60*60*24*7));
    return( "$v weeks, $r" );
}

sub ptime {
    my $msg = shift;
    my $t0 = shift;
    my $out = shift;

    my $timestr= tdiff(time() - $t0);
    my $t1 = time();
    print $out "$msg<font size=\"-1\"> took: $timestr</font>";
    return $t1;
}

sub inflate {
    my $msg = shift;
    my $zf = shift;
    my $tf = shift;
    my $pf = shift;
   
    open(my $hh, ">>$pf"); # or..
    print $hh "<LI><I>$msg:</I>&nbsp;"; 
    my $cmd = "file $zf";
    my $res = qx/$cmd/;
    if ($_debug) {
        print $hh "<BR>$cmd: $res<BR>";
    }
    close $hh; # or..
    open($hh, ">>$pf"); # or..
    if ($res =~ m/gzip compressed data/) {
        $time_t1 = time();
        gunzip $zf => $tf;
        rename $tf => $zf;
        ptime("gunzipped. ", $time_t1, $hh);
    } else {
        print $hh "not gzip compressed.";
    }   
    print $hh "</LI>\n";
    close $hh; # or..
}

sub printIDs {
    my $ids = shift;
    my $of = shift;
    my $out = shift;

    open(my $ofh, ">$of");
    my @INPUT= split("[\n]",$ids);
    my $n=0;
    foreach (@INPUT) {
	$n++;
	$_ =~ s/\x0D$//; # f.ing windows/unix compatibility pb
	$_ =~ s/^[ \t\n]*//; # delete leading/trailing whitespace (spaces, tabs, empty lines)
	$_ =~ s/[ \t\n]*$//;
	if ($_=~/./) { # skip empty lines
# allow whitespaces in query
#	    if ($_ =~ /[ \t]/) {
#               print $out "<b>ERROR: Input format error.</b><br> A list of gene/transcript/protein identifiers requires only one entry per line<br> Please check the line $n, entry= $_\n";
#           }
#        ($_debug) && print "input: $_ <br>";
        if($_debug) {
	    print "$_\n";
	}
        print $ofh "$_\n"
    	}
    }
    close $ofh;
    return $n;
}

sub pprint {
    my $ff = shift;
    my $msg = shift;

    open(my $hh, ">>$ff"); 
    print $hh $msg;
    close($hh);
}
