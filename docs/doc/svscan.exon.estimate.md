<h2> How to estimate fusion gene partner exons</h2>

|title|explanation|
|-----|-----------
|contig catenation| the catenation manner of reference contigs of hgene and tgene|
|hstrand|which strand hgene is on reference
|tstrand|which strand tgene is on reference
|hintron|intron number hgene breakpoint is at
|tintron|intron number tgene breakpoint is at
|hexon|estimated exon number of hgene fused
|texon|estimated exon number of tgene fused

<table>
     <tr>
        <td><h3>contig catenation</h3></td>
        <td><h3>hstrand</h3></td>
        <td><h3>tstrand</h3></td>
        <td><h3>hintron</h3></td>
        <td><h3>tintron</h3></td>
        <td><h3>hexon</h3></td>
        <td><h3>texon</h3></td>
    </tr>
    <tr>
        <td rowspan=4>5'->5'</td>
        <td>+</td>
        <td>+</td>
        <td>hi</td>
        <td>ti</td>
        <td>hi</td>
        <td>ti</td>
    </tr>
    <tr>
        <td>+</td>
        <td>-</td>
        <td>hi</td>
        <td>ti</td>
        <td>hi</td>
        <td>ti+1</td>
    </tr>
    <tr>
        <td>-</td>
        <td>-</td>
        <td>hi</td>
        <td>ti</td>
        <td>hi+1</td>
        <td>ti+1</td>
    </tr>   
    <tr>
        <td>-</td>
        <td>+</td>
        <td>hi</td>
        <td>ti</td>
        <td>hi+1</td>
        <td>ti</td>
    </tr>
    <tr>
        <td rowspan=4>3'->3'</td>
        <td>+</td>
        <td>+</td>
        <td>hi</td>
        <td>ti</td>
        <td>hi+1</td>
        <td>ti+1</td>
    </tr>
    <tr>
        <td>+</td>
        <td>-</td>
        <td>hi</td>
        <td>ti</td>
        <td>hi+1</td>
        <td>ti</td>
    </tr>
    <tr>
        <td>-</td>
        <td>-</td>
        <td>hi</td>
        <td>ti</td>
        <td>hi</td>
        <td>ti</td>
    </tr>  
    <tr>
        <td>-</td>
        <td>+</td>
        <td>hi</td>
        <td>ti</td>
        <td>hi</td>
        <td>ti+1</td>
    </tr>
    <tr>
        <td rowspan=4>5'->3'</td>
        <td>+</td>
        <td>+</td>
        <td>hi</td>
        <td>ti</td>
        <td>hi</td>
        <td>ti+1</td>
    </tr>
    <tr>
        <td>+</td>
        <td>-</td>
        <td>hi</td>
        <td>ti</td>
        <td>hi</td>
        <td>ti</td>
    </tr>
    <tr>
        <td>-</td>
        <td>-</td>
        <td>hi</td>
        <td>ti</td>
        <td>hi+1</td>
        <td>ti</td>
    </tr>
    <tr>
        <td>-</td>
        <td>+</td>
        <td>hi</td>
        <td>ti</td>
        <td>hi+1</td>
        <td>ti+1</td>
    </tr>
    <tr>
        <td rowspan=4>3'->5'</td>
        <td>+</td>
        <td>+</td>
        <td>hi</td>
        <td>ti</td>
        <td>hi+1</td>
        <td>ti</td>
    </tr>
    <tr>
        <td>+</td>
        <td>-</td>
        <td>hi</td>
        <td>ti</td>
        <td>hi+1</td>
        <td>ti+1</td>
    </tr>
    <tr>
        <td>-</td>
        <td>-</td>
        <td>hi</td>
        <td>ti</td>
        <td>hi</td>
        <td>ti+1</td>
    </tr>
    <tr>
        <td>-</td>
        <td>+</td>
        <td>hi</td>
        <td>ti</td>
        <td>hi</td>
        <td>ti</td>
    </tr>
</tbale>


<p> if one partner breakpoint is in exon, then the estimated exon is just that exon and skip the estimation of this partner, but the other partner will undergo estimation if its breakpoint is in intron</p>
