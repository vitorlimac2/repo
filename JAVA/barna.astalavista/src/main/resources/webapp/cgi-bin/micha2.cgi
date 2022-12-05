#!/usr/bin/perl -w

use strict;
use CGI;

print "Content-type: text/html\n\n";

my $cgi = new CGI;


    my $io_handle = $cgi->upload('file_upload');
    my $some_file_path = "../jobs/test.txt";
    print "copying $io_handle to $some_file_path\n";
    open ( my $out_file, ">$some_file_path" ) or die $!;
    print "while";
    while ( my $bytesread = $io_handle->read(my $buffer,1024) ) {
        print $out_file $buffer;
        print $buffer;
    }
    close ($out_file);

