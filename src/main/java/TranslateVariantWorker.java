package main.java;


import org.apache.commons.lang3.StringUtils;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static main.java.DbCreator.translateDNA;


public class TranslateVariantWorker implements  Runnable{


    private JVariant jv = new JVariant();
    private boolean verbose = false;

    public static HashMap<String,JTranscript> gMap = new HashMap<>();

    public static String prefix4variantID = "VAR";

    public BufferedWriter logWriter = null;

    private int nvar = 0;

    public TranslateVariantWorker(JVariant jVariant, int nv, boolean vb,BufferedWriter logWriter) {
        this.jv = jVariant;
        this.nvar = nv;
        this.verbose = vb;
        this.logWriter = logWriter;
    }


    @Override
    public void run() {
        String ID = jv.ID;
        String transcript = jv.transcriptID;
        int start = jv.start;
        int end = jv.end;
        String ref = jv.ref;
        String obs = jv.obs;
        String cchange = jv.cchange;
        String pchange = jv.pchange;
        String chr = jv.chr;

        if(verbose){
            try {
                logWriter.write("Process: "+jv.ID+" with "+jv.transcriptID+", "+jv.start+", "+jv.end+", "+jv.ref+", "+jv.obs+", "+jv.cchange+", "+jv.pchange+"\n");
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        if(!gMap.containsKey(jv.transcriptID)){
            try {
                logWriter.write("WARNING: cannot find mRNA sequence for "+jv.transcriptID+" in the fasta file!\n");
            } catch (IOException e) {
                e.printStackTrace();
            }
            jv.valid = false;
            return;
        }
        if(!gMap.containsKey(jv.transcriptID)){
            try {
                logWriter.write("WARNING: cannot find annotation for "+jv.transcriptID+" in the gene file or cannot infer the transcription start site\n");
            } catch (IOException e) {
                e.printStackTrace();
            }
            jv.valid = false;
            return;
        }

        if(jv.end > gMap.get(jv.transcriptID).mRNA_seq.length()){
            try {
                logWriter.write("ERROR: transcript end ("+gMap.get(jv.transcriptID).end+") for "+jv.transcriptID+" is longer than transcript length "+gMap.get(jv.transcriptID).mRNA_seq.length()+", skipping this transcript\n");
            } catch (IOException e) {
                e.printStackTrace();
            }
            jv.valid = false;
            return;
        }

        // not used
        if(false){
            String utr5 = gMap.get(transcript).mRNA_seq.substring(0,gMap.get(transcript).start-1);
            String utr3 = gMap.get(transcript).mRNA_seq.substring(gMap.get(transcript).end);
            String mrna1;
            String mrna2;
            mrna1 = gMap.get(transcript).mRNA_seq;

            String dna_seq = gMap.get(jv.transcriptID).mRNA_seq.substring(gMap.get(jv.transcriptID).start-1, gMap.get(jv.transcriptID).end - gMap.get(jv.transcriptID).start + 1);
            String dna[] = dna_seq.split("");
            ArrayList<String> dnaList = new ArrayList<>();
            for(int i=0;i<dna.length;i++){
                dnaList.add(dna[i]);
            }
            String protein1 = null;
            String protein2;
            String warning = "";

            String cds = dna_seq;

            if(jv.end > dna.length){
                try {
                    logWriter.write("ERROR in "+jv.ID+": end position of variant ("+jv.end+") in "+jv.transcriptID+" is longer than coding portion length "+dna.length+", skipping this transcript\n");
                } catch (IOException e) {
                    e.printStackTrace();
                }
                jv.valid = false;
                return;
            }

            try {
                protein1 = translateDNA(dna_seq,logWriter);
            } catch (IOException e) {
                e.printStackTrace();
            }
            Pattern pPro = Pattern.compile("\\*\\w");
            Matcher matcher = pPro.matcher(protein1);
            if(!matcher.find()){
                warning = "(WARNING: Potential FASTA sequence error!!!)";
            }

            // we need to confirm here
            if(jv.start == jv.end && jv.ref.isEmpty() && !jv.obs.isEmpty()){
                // this is an insertion
                dnaList.add(jv.start,jv.obs);
            }else{
                for(int i=jv.start;i<(jv.end-jv.start+1);i++){
                    dnaList.remove(i);
                }
                dnaList.add(jv.start,jv.obs);
            }

            String dna_seq_variant = StringUtils.join(dnaList,"");

            mrna2 = utr5 + dna_seq_variant + utr3;


        }else{

            String dna_seq = gMap.get(transcript).mRNA_seq.substring(gMap.get(transcript).start-1, gMap.get(transcript).end);


            String dna_seq_utr = gMap.get(transcript).mRNA_seq.substring(gMap.get(transcript).start-1);
            String dna[] = dna_seq_utr.split("");
            ArrayList<String> dnaList = new ArrayList<>();
            for(int i=0;i<dna.length;i++){
                dnaList.add(dna[i]);
            }
            String protein1 = null;
            String protein2 = null;

            String warning = "";
            int inframe = 0;

            if(jv.end > dna.length){
                try {
                    logWriter.write("WARNING in "+ID+": end position of variant ("+end+") in "+transcript+" is longer than coding portion length "+dna.length+", possibly due to a deletion on a stop codon\n");
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }

            if(verbose){
                try {
                    logWriter.write(gMap.get(transcript).mRNA_seq+"\t"+(gMap.get(transcript).start-1)+","+(gMap.get(transcript).end - gMap.get(transcript).start + 1)+"\n");
                } catch (IOException e) {
                    e.printStackTrace();
                }
                try {
                    logWriter.write("NOTICE: wild-type DNA sequence is "+dna_seq+"\n");
                } catch (IOException e) {
                    e.printStackTrace();
                }
                try {
                    logWriter.write("Include UTR: "+dna_seq_utr+"\n");
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }

            try {
                protein1 = translateDNA(dna_seq,logWriter);
            } catch (IOException e) {
                e.printStackTrace();
            }

            // System.out.println(protein1);

            Pattern pStop = Pattern.compile("\\*.+");
            Matcher pStopMatcher = pStop.matcher(protein1);

            if(pStopMatcher.find()){

                warning = "(WARNING: Potential FASTA sequence error!!!)";
                jv.valid = false;
                return;
            }

            if(start == end && ref.isEmpty() && !obs.isEmpty()){
                // this is an insertion (ref is treated as 0 when we read the EVF file previously)
                dnaList.add(start,obs);
                start = start + 1;
                end = end + 1;
                if(obs.length() % 3 == 0){
                    inframe = inframe + 1;
                }

            }else { // this is a substitution
                if(obs.isEmpty()){ // deletion
                    // AAA -> ""
                    // fromIndex, inclusive, and toIndex, exclusive.
                    dnaList.subList(start-1,end).clear();
                    if((end - start + 1) % 3 == 0){
                        inframe = inframe + 1;
                    }

                }else {
                    // absolutely important here to split obs into individual elements in an array!!!
                    dnaList.subList(start-1,end).clear();
                    // Inserts the specified element at the specified position in this
                    // list. Shifts the element currently at that position (if any) and
                    // any subsequent elements to the right (adds one to their indices).
                    // index index at which the specified element is to be inserted
                    dnaList.add(start-1,obs);


                    if((end - start + 1 - obs.length() % 3 == 0) ){
                        inframe = inframe + 1;
                    }


                }
            }

            int aastart = (int) ((start-1)/3) + 1 ;
            int aaend1 = (int)((end-1)/3)+1;
            int aaend2 = (int)((start+obs.length()-1-1)/3) + 1;
            String aachange = ""; // assign default value, since synonymous variants do not have aachange
            String dna_seq_variant = StringUtils.join(dnaList,"");

            // protein2 should be truncated so that anything after stop codon be deleted
            try {
                protein2 = translateDNA(dna_seq_variant,logWriter);
            } catch (IOException e) {
                e.printStackTrace();
            }

            protein2 = protein2.replaceAll("\\*.+","*");

            if(verbose){

            }

            String function;

            String protein2_tmp;
            if(aaend2 <= protein2.length()) {
                protein2_tmp = protein2.substring(aastart - 1, aaend2);
            }else{
                protein2_tmp = protein2.substring(aastart - 1);
            }

            Pattern pStopOnly = Pattern.compile("\\*");
            Matcher pStopOnlyMatcher = pStopOnly.matcher(protein2_tmp);
            if(pStopOnlyMatcher.find()){
                function = "immediate-stopgain";
            }else if(end >= (gMap.get(transcript).end-2)){
                function = "immediate-stoploss";
            }else {
                if(protein1.equalsIgnoreCase(protein2)){
                    function = "silent";
                    if(verbose){
                        try {
                            logWriter.write("NOTICE: silent variant!\n");
                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                        jv.valid = false;
                        return;
                    }
                }else{
                    function = "protein-altering"; // change nonsynonymous to protein-altering to reduce confusion among users
                }
            }


            // The goal is to find the exact change of all amino acids (not just the first one), until a stop codon is found
            // pos1: first position with difference aa1:peptide with difference

            int pos1 = 0;
            int pos2 = 0;
            String aa1 = "";
            String aa2 = "";
            int len_max = protein1.length();
            if(len_max < protein2.length()){
                len_max = protein2.length();
            }

            String pro1List[] = protein1.split("");
            String pro2List[] = protein2.split("");

            for(int i=0;i<len_max;i++){
                if(pos1>=1){

                    if(i < protein1.length()){
                        aa1 = aa1 + pro1List[i];
                    }
                    if(i < protein2.length()){
                        aa2 = aa2 + pro2List[i];
                    }
                }else{
                    if(i > protein1.length()){
                        continue;
                    }

                    if(i > protein2.length()){
                        continue;
                    }

                    if(i == protein1.length()){
                        pos1 = i+1;
                        aa1 = "";
                        aa2 = pro2List[i];
                    }else if(i == protein2.length()){
                        pos1 = i+1;
                        aa1 = pro1List[i];
                        aa2 = "";
                    }else if(!pro1List[i].equalsIgnoreCase(pro2List[i])){
                        pos1 = i+1;
                        aa1 = pro1List[i];
                        aa2 = pro2List[i];
                    }
                }

            }

            // maximum of pos2 is the end of the wildtype protein
            if(pos2==0){
                pos2 = protein1.length();
            }


            String aa1d[] = aa1.split("");
            String aa2d[] = aa2.split("");
            ArrayList<String> aa1List = new ArrayList<>();
            ArrayList<String> aa2List = new ArrayList<>();
            for(int i=0;i<aa1d.length;i++){
                aa1List.add(aa1d[i]);
            }

            for(int i=0;i<aa2d.length;i++){
                aa2List.add(aa2d[i]);
            }


            // sometimes there is a nonframeshift mutation, and I previously calculate the protein until stop
            // codon, which is inconvenient.
            // the following code addresses this problem.
            while(true){  // sometimes aa1 and aa2 may have identical tails

                if(aa1List.size()>=1 && aa2List.size()>=1) {
                    String curAA1 = StringUtils.join(aa1List, "");
                    String lastAA1 = aa1List.get(aa1List.size() - 1);
                    String lastAA2 = aa2List.get(aa2List.size() - 1);

                    // System.out.println("aa1list size:" + aa1List.size() + ", aa2list size:" + aa2List.size() + "\t" + lastAA1 + "\t" + lastAA2);

                    if (pos1 >= 1 && !curAA1.equalsIgnoreCase("*") && lastAA1.equalsIgnoreCase(lastAA2)) {
                        aa1List.remove(aa1List.size() - 1);
                        aa2List.remove(aa2List.size() - 1);
                        pos2 = pos2 - 1;
                    } else {
                        break;
                    }
                }else{
                    break;
                }
            }

            aa1 = StringUtils.join(aa1List,"");
            aa2 = StringUtils.join(aa2List,"");

            // protein position
            String position = "-";

            if(pos1 == 0 ){
                aachange = " (no amino acid change)";

                Pattern fsPattern = Pattern.compile("([\\w\\*]+?)(\\d+)fs");
                Matcher fsMatcher = fsPattern.matcher(pchange);
                if(pchange.isEmpty()){
                    // it is actually not that simple to handle this situation of silent change. We have to calculate it from c change instead
                    Pattern numPattern = Pattern.compile("(\\d+)");
                    Matcher numMatcher = numPattern.matcher(cchange);
                    if(numMatcher.find()){
                        int pchange_pos = (int)((Integer.valueOf(numMatcher.group(1)) -1)/3);
                        if(pchange_pos >= protein1.length()){ // something is wrong here, the wildtype sequence is not complete possibly

                            jv.valid = false;
                            return;
                        }
                        String pchange_aa = protein1.substring(pchange_pos-1,pchange_pos);
                        pchange = pchange_aa+pchange_pos+pchange_aa;
                        position = String.valueOf(pchange_pos);
                    }
                }else if(fsMatcher.find()){
                    pchange = fsMatcher.group(1)+fsMatcher.group(2)+fsMatcher.group(1);
                    position = fsMatcher.group(2);
                }

                function = "silent";

            }else if(pos1 == pos2){ // a simple non-synonymous change
                // a simple non-synonymous change, or an inframe insertion if length of aa1!=aa2

                aachange = " (position "+pos1+" changed from "+aa1+" to "+aa2+")";
                pchange = aa1+pos1+aa2;
                position = String.valueOf(pos1);

                if(aa1.equalsIgnoreCase("*")){
                    function = "immediate-stoploss";
                }
                if(aa2.equalsIgnoreCase("*")){
                    function = "immediate-stopgain";
                }

                if(aa2.isEmpty()){
                    pchange = pchange + "del";
                }else if(aa1.length() != aa2.length()){
                    pchange = aa1 + pos1 + "delins" + aa2;
                }

            }else if(pos1 > pos2){ // when there is an insertion, pos2 will be 1 less than pos1 (example: 1	1647893 1647893 -       TTTCTT)
                aachange = " (position "+pos2+"-"+pos1+" has insertion "+aa2+")";
                position = pos2+"-"+pos1;

                if(pos1<=protein1.length()){
                    pchange = protein1.substring(pos2-1,pos2)+pos2+"_"+protein1.substring(pos1-1,pos1)+pos1+"ins"+aa2;
                }else{
                    pchange = protein1.substring(pos2-1,pos2)+pos2+"_"+pos1+"ins"+aa2;
                }


            }else{

                position = pos1+"-"+pos2;

                if(inframe>0){ // inframe substitution
                    aachange = " (position "+pos1+"-"+pos2+" changed from "+aa1+" to "+aa2+")";
                    pchange = aa1.substring(0,1)+pos1+"_"+aa1.substring(aa1.length()-1)+pos2+"delins"+aa2;
                    if(aa2.isEmpty()){
                        // when aa2 is empty, it should be just deletion, not delins
                        pchange.replaceAll("ins$","");
                    }
                }else {
                    aachange = " (position "+pos1+"-"+pos2+" changed from "+aa1+" to "+aa2+")";
                    if(aa2.length()==0) {
                        pchange = aa1.substring(0, 1) + pos1 + "fs*" + aa2.length();
                    }else{
                        pchange = aa1.substring(0, 1) + pos1 + aa2.substring(0, 1) + "fs*" + aa2.length();
                    }
                }
            }

            protein2 = protein2.replaceAll("\\*$","");


            if(function.equalsIgnoreCase("silent")){
                try {
                    logWriter.write("Warning: silent => "+ transcript+"\tc."+cchange+"\tp."+pchange+"\t"+aachange+"\n");
                } catch (IOException e) {
                    e.printStackTrace();
                }
                jv.valid = false;
                return;
            }

            jv.variant_function = function;

            // nvar++;
            // ID:
            String new_variant_protein_ID = prefix4variantID+"|"+transcript+"|"+nvar;


            if(aa2.isEmpty()){
                if(jv.variant_function.equalsIgnoreCase("immediate-stopgain")){
                    aa2 = "*";
                }else{
                    aa2 = "-";
                }
            }

            StringBuilder vBuilder = new StringBuilder();
            vBuilder.append(new_variant_protein_ID).append("\t")
                    .append(chr).append("\t")
                    .append(jv.genome_start).append("\t")
                    .append(jv.genome_end).append("\t")
                    .append(jv.genome_ref).append("\t")
                    .append(jv.genome_var).append("\t")
                    .append(jv.variant_type).append("\t")
                    .append(jv.variant_function).append("\t")
                    .append(jv.geneID).append("\t")
                    .append(jv.transcriptID).append("\t")
                    .append(jv.transcriptID).append("\t") // protein id
                    .append(jv.cchange).append("\t")
                    .append(jv.pchange).append("\t")
                    .append(aa1).append("\t")
                    .append(position).append("\t")
                    .append(aa2);
            jv.valid = true;
            jv.ID = new_variant_protein_ID;
            jv.line = vBuilder.toString();
            jv.protein_seq = protein2;
            jv.protein_header = transcript+" c."+cchange+" p."+pchange+" "+function+" "+aachange;
            // vWriter.write(vBuilder.toString());

        }

    }
}
