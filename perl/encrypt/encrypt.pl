#!/usr/bin/perl
use strict;
use warnings;
use Crypt::RC4::XS;
use MIME::Base64;

open my $fh, '<', 'snpinfo_171224_ori.txt' or die "Can't open file $!";
my $plaintext = do { local $/; <$fh> };

my $key       = "pzweuj";
my $encrypted = RC4($key, $plaintext);

my $encoded = encode_base64($encrypted);
print $encoded
