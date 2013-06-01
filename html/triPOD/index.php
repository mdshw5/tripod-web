
<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd"> 
 
<html> 
<!--
########################################
#  Copyright (c)  2012 Joseph Baugher, Matthew Shirley and Jonathan Pevsner.
#                 All Rights Reserved.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS "AS IS" AND ANY 
#  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
#  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNERS BE
#  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
#  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
#  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
#  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
#  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
#  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
#  POSSIBILITY OF SUCH DAMAGE. THIS SOFTWARE IS FREE FOR PERSONAL OR ACADEMIC
#  USE. THE SOFTWARE MAY NOT BE USED COMMERCIALLY WITHOUT THE EXPRESS, WRITTEN
#  PERMISSION OF THE COPYRIGHT HOLDERS. ALL ACTIONS OR PROCEEDINGS RELATED TO 
#  THIS SOFTWARE SHALL BE TRIED EXCLUSIVELY IN THE STATE AND FEDERAL COURTS 
#  LOCATED IN THE COUNTY OF BALTIMORE CITY, MARYLAND. USE OF THIS SOFTWARE
#  IMPLIES ACCEPTANCE OF THIS CONTRACT.
########################################
--> 
<head> 
<link rel="stylesheet" type="text/css" href="triPOD.css">
<title>triPOD at the Pevsnerlab -- Upload Page</title> 
<meta name="author" content="Matt Shirley" > 
<meta name="description" content="Upload page for the triPOD web tool"> 
<meta name="keywords" content="chromosomal,abnormalities,mosaic,SNP,microarray,parent,origin,trio,BAF"> 
</head> 
 
<body> 
<script type="text/javascript">

   var _gaq = _gaq || [];
_gaq.push(['_setAccount', 'UA-29299688-1']);
_gaq.push(['_trackPageview']);

(function() {
  var ga = document.createElement('script'); ga.type = 'text/javascript'; ga.async = true;
  ga.src = ('https:' == document.location.protocol ? 'https://ssl' : 'http://www') + '.google-analytics.com/ga.js';
  var s = document.getElementsByTagName('script')[0]; s.parentNode.insertBefore(ga, s);
})();

</script>

<?php include 'header.php'; ?>

<form enctype="multipart/form-data" action="/cgi-bin/triPOD/handler.py" method="post"> 
 
<table width="100%" border="0" cellpadding="6" cellspacing="3"> 
 
<tr class="rowcommon"> 
<td class="cell1 cell1width">Step 1</td> 
<td class="cell2 cell2width">Select a file to upload</td> 
<td class="cell2"><input type="file" name="file" size="25"><br><input type="checkbox" name="sampledata" value="sampledata.txt">Use sample data</td> 
</tr> 

<tr class="rowcommon"> 
<td class="cell1 cell1width">Step 2</td> 
<td class="cell2 cell2width">Statistical significance threshold for p-value</td> 
<td class="cell2"><input name="alpha"  type="text" size="10" value="0.1"></td> 
</tr> 

<tr class="rowcommon"> 
<td class="cell1 cell1width">Step 3</td> 
<td class="cell2 cell2width">Reference genome build</td> 
<td class="cell2"> 
<select name="build"> 
<option value="hg16_centromeres.txt">hg16</option> 
<option value="hg17_centromeres.txt">hg17</option> 
<option value="hg18_centromeres.txt" selected>hg18</option> 
<option value="hg19_centromeres.txt">hg19(b37)</option> 
</select> 
</td> 
</tr> 

<tr class="rowcommon"> 
<td class="cell1 cell1width">Step 4</td> 
<td class="cell2 cell2width">Gender of sample</td> 
<td class="cell2"> 
<select name="gender"> 
<option value="M">Male</option> 
<option value="F">Female</option> 
<option value="NA" selected>NA</option> 
</select> 
</td> 
</tr> 
 
<tr class="rowcommon"> 
<td class="cell1 cell1width">Step 5</td> 
<td class="cell2 cell2width">Submit query. Note that large files will take some time to upload.</td> 
<td class="cell2"><input type="submit"></td> 
</tr>  

<tr class="rowcommon">
<td class="cell1 cell1width">(optional) Detection Methods</td> 
<td class="cell2 cell2width">Changing the detection methods used will alter your results.</td>
<td class="cell2">
<a id="displayText" href="javascript:toggle();">Show detection methods</a>

</td>

</table>

<div id="toggle" style="display: none">
<table width="100%" border="0" cellpadding="6" cellspacing="3"> 

<tr class="rowcommon"> 
<td class="cell1 cell1width">POD</td> 
<td class="cell2 cell2width">Parent of origin detection method</td> 
<td class="cell2"> 
<select name="pod"> 
<option value="pod" selected>Yes</option> 
<option value="nopod">No</option> 
</select> 
</td> 
</tr> 

<tr class="rowcommon"> 
<td class="cell1 cell1width">PODhd</td> 
<td class="cell2 cell2width">Homozygous deletion detection method</td> 
<td class="cell2"> 
<select name="podhd"> 
<option value="hd">Yes</option> 
<option value="nohd" selected>No</option> 
</select> 
</td> 
</tr> 

<tr class="rowcommon"> 
<td class="cell1 cell1width">PODmi1</td> 
<td class="cell2 cell2width">Single Mendelian error detection method</td> 
<td class="cell2"> 
<select name="podmi1"> 
<option value="mi1">Yes</option> 
<option value="nomi1" selected>No</option> 
</select> 
</td> 
</tr> 

<tr class="rowcommon"> 
<td class="cell1 cell1width">PODcr</td> 
<td class="cell2 cell2width">Outlier B allele frequency detection method</td> 
<td class="cell2"> 
<select name="podcr"> 
<option value="podcr">Yes</option> 
<option value="nopodcr" selected>No</option> 
</select> 
</td> 
</tr> 

</table> 
</div>
 
</form> 
<a href="triPOD_simulation_trios.zip">Download simulated benchmarking data (zip)</a><br>
<a href="triPOD_simulation_trios.tar.bz2">Download simulated benchmarking data (bzip2)</a><p>
   <a style="color:gray">This web application is currently under active development. If unexpected errors are encountered, please try using a browser other than Internet Explorer, or a new data file.</a>
</body> 
 
<script language="javascript"> 
function toggle() {
	var ele = document.getElementById("toggle");
	var text = document.getElementById("displayText");
	if(ele.style.display == "block") {
    		ele.style.display = "none";
		text.innerHTML = "Show detection methods";
  	}
	else {
		ele.style.display = "block";
		text.innerHTML = "Hide detection methods";
	}
} 
</script>
 
</html> 
