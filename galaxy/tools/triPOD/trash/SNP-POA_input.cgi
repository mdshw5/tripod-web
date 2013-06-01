#!/usr/bin/perl -w
use CGI qw/:standard/;
use CGI::Carp qw/fatalsToBrowser/; # this will print all of the errors to the browser for debugging
use Digest::MD5  qw(md5_hex);
use strict;
## configuration
# web server root
my $webroot = "/var/www";
# path to the html root directory to store files
my $htmlroot = $webroot."/html/SNP-POA";
# path to the directory to store results
my $resultsroot = $webroot."/html/SNP-POA/results";
# path to the cgi-bin directory to store this program
my $cgiroot = $webroot."/cgi-bin/SNP-POA";
# relative path to SNP-POA results directory
my $webresults = "../../SNP-POA/results";
# path to perl binaries
my $perlroot = "/usr/bin/perl";

# CGI form input
my $time = time();
my $query = CGI->new;
my $cores = $query->param('cores');
my $bin = $query->param('sizeRegion');
my $prob = $query->param('significanceThreshold');
my $hetSD = $query->param('numHetStDevs');
my $homSD = $query->param('numHomzygStDevs');
my $errNC = $query->param('ErrorAsPercentNC');
my $err = $query->param('ErrorRate');
my $nc = $query->param('NoCallThreshold');
my $gender = $query->param('gender');
my $filename = $query->param('file');

# In case there is a problem with the file upload we want to print an error message to the browser about what happened.
   if (!$filename && $query->cgi_error) {
      print $query->header(-status=>$query->cgi_error);
      exit 0;
   }

# Grab the full path to the temporary file that CGI module creates, so we don't have to store any extra uploaded files
my $tmpfilename = $query->tmpFileName($filename);

# Generate a unique hash based on the date and time, so that we don't overwrite other results
my $md5_hash = md5_hex(localtime());
my $resultsdir = join("/", $resultsroot, $md5_hash); 

# SNP-POA system call configuration 
my $SNPPOA = join(" ", "/var/www/cgi-bin/SNP-POA/SNP-POA.pl",$tmpfilename,
	"--out",$resultsdir,
	"--cores",$cores,
	 "--bin",$bin,
	"--prob",$prob,
	"--hetSD",$hetSD,
	"--homSD",$homSD,
	"--err",$errNC,
	"--erate",$err,
	"--nc",$nc,
	"--gender",$gender,
	 "1>> /var/log/SNP-POA/log");

my @error;
if (defined($tmpfilename)) { @error = qx($SNPPOA); }
	else { 
		print $query->header(), $query->start_html(), "Please select a file to upload", $query->end_html(); 
		die();
}
# Catch errors from executing $SNPPOA
if (defined(@error)) { 
		print $query->header(), $query->start_html(), @error, $query->end_html();
                die();
}

# We'll make a .zip archive of the results directory
my $archivename = "SNP-POA.zip";
if (-d $resultsdir) {
  chdir $resultsdir;
  system(join(" ","zip","-r",$archivename,"*","1> /dev/null"));
}
else {
	print $query->header(), $query->start_html(), "Could not create archive file", $query->end_html();
	die();
}

# Start gathering things we need to make the results HTML
open(REGIONS,"<".join("/",$resultsroot,$md5_hash,"regionFile.txt"));
  my @regions = <REGIONS>;
close(REGIONS);
my $archiveurl = join("/", $webresults,$md5_hash,$archivename);
my $imagepath = join("/", $resultsroot, $md5_hash, "Plots/");
my $webimagepath = join("/", $webresults, $md5_hash, "Plots/");
my @images = <$imagepath*>; # Get a listing of all image files in the results directory

## begin HTML output
print 
  $query->header(),
  $query->start_html('SNP-POA Results'),
  "Finished in: ", time() - $time, " seconds.", "<p>",
  "<a>Archive (zip) containing results: </a>",
  "<a href=$archiveurl>$archivename</a>", "<p>";
print "<table align=\"left\" border=\"1\">";
# print out a formatted table from the regions file
foreach (@regions) {
  print "<tr>";
  my @line = split(/[\t]/,$_);
    foreach(@line) {
      print "<td>$_</td>";
    }
  print "</tr>";
}
#
# print images in the same table
foreach (@images) {
  my @imagename = split(/[\/]/,$_);
  print 
    "<tr><td>
    <img src=",join("/", $webimagepath, $imagename[-1])," width=\"640\" align=\"center\">
    </td></tr>";
}
print "</table>";
print
  "<p>",
  $query->end_html();
exit;
