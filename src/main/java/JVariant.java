package main.java;

import java.util.ArrayList;

/**
 * sava a variant record
 */
public class JVariant {

    /**
     * Variant ID
     */
    public String ID;

    public int start;
    public int end;
    public String ref;
    public String obs;
    public String transcriptID;
    public String geneID;
    public String cchange;
    public String pchange;
    public String proteinID;
    public String chr;

    public int genome_start;
    public int genome_end;
    public String variant_type;

    public String variant_function;

    public String genome_ref = "";
    public String genome_var = "";

    public String protein_seq;

    /**
     * The information will be exported in fasta header position
     */
    public String protein_header = "";
    /**
     * This is used to save the information of output string.
     */
    public String line = "";
    public static String lineHead = "";
    public ArrayList<String> samples = new ArrayList<>();

    public boolean valid = false;

}
