#!/usr/bin/perl
use strict;
use warnings;
use Crypt::RC4::XS;
use MIME::Base64;

open my $fh, '<', 'inputsnp' or die "Can't open file $!";
my $plaintext = do { local $/; <$fh> };

my $decoded = decode_base64($plaintext);
my $key       = "pzweuj";
my $decrypted = RC4($key, $decoded);

# print $decrypted

open my $snpinputdecry, '<', \$decrypted;
while (<$snpinputdecry>) {
  print "$_";
}