package main.java;

public class JTranscript {

    public int start;
    public int end;
    public String transcriptID;

    /**
     * can be used to distinguish whether it's coding gene or not.
     */
    public int cds_start = 0;
    public int cds_end = 0;

    /**
     * DNA sequence from hg19_refGeneMrna.fa
     */
    public String mRNA_seq = "";

    /**
     *
     * This variable is used to save the information for confirming the transcript in hg19_refGeneMrna.fa and
     * hg19_refGene.txt. Sometimes, there are multiple rows from the same transcript. In this situation, we need to
     * select the true transcript based on the following information from hg19_refGeneMrna.fa
     * Comment: this sequence (leftmost exon at chr1:222001007) is generated ...
     * The format is chr1:222001007
     */
    public String infor = "";
}
