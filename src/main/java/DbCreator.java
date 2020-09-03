package main.java;

import net.sf.jfasta.FASTAElement;
import net.sf.jfasta.FASTAFileReader;
import net.sf.jfasta.impl.FASTAElementIterator;
import net.sf.jfasta.impl.FASTAFileReaderImpl;
import org.apache.commons.cli.*;
import org.apache.commons.lang3.StringUtils;

import java.io.*;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Build a customized protein database
 */
public class DbCreator {

    public static HashMap<String,String> codons;


    /**
     * The prefix of variant protein ID, default is VAR_
     */
    public String prefix4variantID = "VAR";


    public String output_var_table_suffix = "-varInfo.txt";

    public String output_db_table_suffix = "-var.fasta";

    public String output_stat_table_suffix = "-varStat.txt";

    /**
     * mRNA annotation data and mRNA sequences
     */
    private HashMap<String,JTranscript> gMap = new HashMap<>();


    public String outdir = "./";

    /**
     * Whether or not to add reference protein sequences separately
     */
    public String outputRefProteinDB = "";


    /**
     * Construct a new DbCreator object
     * @param mRNA_db *_refGeneMrna.fa
     * @param gene_anno_file *__refGene.txt
     */
    public DbCreator(String mRNA_db, String gene_anno_file, BufferedWriter logWriter,String out_dir){
        this.outdir = out_dir;
        try {
            this.readAnnotationData(mRNA_db,gene_anno_file,logWriter);
        } catch (IOException e) {
            e.printStackTrace();
        }
        codons = getCodon();
    }


    public static void main(String[] args) throws ParseException, IOException {
        Options options = new Options();

        // input file
        // sample somatic germline rna
        // s1 s.txt g.txt r.txt
        // s2 s.txt g.txt r.txt
        // ...

        // output file
        // one sample one column,
        // s1 s2 ...
        // 0,1,1 1,1,1

        options.addOption("i", true,"Annovar annotation result file. Multiple files are separated by ','.");
        options.addOption("f", true,"A file which includes multiple samples. This parameter is used to build a customized database for ");
        options.addOption("r", true,"Gene annotation data");
        options.addOption("d", true,"mRNA fasta database");
        options.addOption("o", true,"Output folder");
        options.addOption("p1", true,"The prefix of variant protein ID, default is VAR_");
        options.addOption("p2", true,"The prefix of final output files, default is merge");
        options.addOption("t",false,"Whether or not to add reference protein sequences to the output database file");
        options.addOption("ref", true, "Output reference protein database file");
        options.addOption("v", false, "Verbose");
        options.addOption("h", false, "Help");


        CommandLineParser parser = new DefaultParser();
        CommandLine cmd = parser.parse(options, args);
        if (cmd.hasOption("h") || cmd.hasOption("help") || args.length == 0) {
            HelpFormatter f = new HelpFormatter();
            System.out.println("java -Xmx2G -jar customprodbj.jar");
            // System.out.println(version);
            f.printHelp("Options", options);
            System.exit(0);
        }

        String outdir = "./";
        if(cmd.hasOption("o")){
            outdir = cmd.getOptionValue("o");
            File OD = new File(outdir);
            if(!OD.isDirectory()){
                OD.mkdirs();
            }
        }

        DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
        Date date = new Date();

        String logfile = outdir+"/log.txt";
        System.out.println("log file:"+logfile);

        BufferedWriter logWriter = new BufferedWriter(new FileWriter(new File(logfile)));
        logWriter.write(dateFormat.format(date)+"\n");
        logWriter.write("Command line:"+StringUtils.join(args," ")+"\n");


        String geneAnno = cmd.getOptionValue("r");
        String genefa = cmd.getOptionValue("d");

        DbCreator createDB4Annovar = new DbCreator(genefa,geneAnno,logWriter,outdir);

        // createDB4Annovar.outdir = outdir;


        boolean includeRefProteins = false;
        if(cmd.hasOption("t")){
            includeRefProteins = true;
        }

        if(cmd.hasOption("p1")){
            createDB4Annovar.prefix4variantID = cmd.getOptionValue("p1");
        }

        boolean verbose = false;
        if(cmd.hasOption("v")){
            verbose = true;
        }

        if(cmd.hasOption("f")) {
            // input is a file which includes multiple samples or multiple variant annotation files for the same sample
            // This is used to build a customized protein database for a TMT/iTRAQ experiment which usually contains
            // multiple samples. This can also be used to build a customized protein database for a sample which contains
            // multiple variant results, for example, a variant calling result from WES data, a variant calling result
            // from RNA-Seq data.
            String input_list = cmd.getOptionValue("f");
            File IFL =  new File(input_list);
            BufferedReader iReader = new BufferedReader(new FileReader(IFL));
            String headLine = iReader.readLine().trim();
            String h[] = headLine.split("\\s+");
            HashMap<String, Integer> headerMap = new HashMap<>();
            ArrayList<String> hList = new ArrayList<>();
            for (int i = 0; i < h.length; i++) {
                headerMap.put(h[i], i);
                hList.add(h[i]);
            }

            // read sample list file
            // please note that a file may not exist
            String line;
            HashMap<String, HashMap<String, ArrayList<String>>> resfileMap = new HashMap<>();
            while ((line = iReader.readLine()) != null) {
                String d[] = line.trim().split("\\s+");
                String sample_name = d[0];
                System.out.println("process file:"+sample_name);
                HashMap<String, ArrayList<String>> smap = new HashMap<>();
                for (int i = 1; i < d.length; i++) {
                    System.out.println(" => "+h[i]+" => "+ d[i]);
                    String prefix = sample_name + "-" + h[i];
                    ArrayList<String> rfile = createDB4Annovar.run(d[i], prefix, false,true, true,logWriter);
                    // h[i] is sample name
                    smap.put(h[i], rfile);

                }
                resfileMap.put(sample_name, smap);
            }

            iReader.close();

            // merge
            System.out.println("Merge ...");
            String oprefix = IFL.getName().replaceAll(".txt","");
            createDB4Annovar.merge(resfileMap, hList, oprefix, includeRefProteins,logWriter);
        }else{

            String annovarRes = cmd.getOptionValue("i");
            // we don't need to merge
            String oprefix;
            if(cmd.hasOption("p2")){
                oprefix = cmd.getOptionValue("p2");
            }else{
                oprefix = "merge";
            }
            createDB4Annovar.run(annovarRes,oprefix, includeRefProteins, verbose,true,logWriter);


        }

        if(cmd.hasOption("ref")){
            createDB4Annovar.outputRefProteinDB = cmd.getOptionValue("ref");
        }else{
            createDB4Annovar.outputRefProteinDB = "";
        }

        if(!createDB4Annovar.outputRefProteinDB.isEmpty()){
            System.out.println("Reference protein db:"+createDB4Annovar.outputRefProteinDB);
            createDB4Annovar.addReferenceProteins(createDB4Annovar.gMap, createDB4Annovar.outputRefProteinDB,logWriter);
        }else{
            System.out.println(createDB4Annovar.outputRefProteinDB);
        }

        logWriter.close();

    }


    public void merge(HashMap<String,HashMap<String,ArrayList<String>>> resfileMap, ArrayList<String> h, String prefix, boolean addRef, BufferedWriter logWriter) throws IOException {

        HashMap<String, JVariant> mergeVarMap = new HashMap<>();

        // varID => sample+head
        HashMap<String,HashSet<String>> varID2sample = new HashMap<>();

        // Get all varIDs across all samples
        ArrayList<String> sampleNames = new ArrayList<>();
        for(String sample_name : resfileMap.keySet()){
            sampleNames.add(sample_name);
            for(int i=1;i<h.size();i++) {
                String fileClass = h.get(i);
                String annofile = resfileMap.get(sample_name).get(fileClass).get(0);
                String varfa = resfileMap.get(sample_name).get(fileClass).get(1);

                logWriter.write("Read " + annofile+"\n");
                logWriter.write("Read " + varfa + "\n");
                HashMap<String, JVariant> varMap = readData(annofile, varfa,logWriter);

                // varID => String varID = d[hMap.get("Chr")]+"_"+d[hMap.get("Start")]+"_"+d[hMap.get("End")]+"_"+d[hMap.get("Ref")]+"_"+d[hMap.get("Alt")]+"_"+d[hMap.get("mRNA")];
                String sample_head = sample_name + " " + fileClass;
                for (String varID : varMap.keySet()) {
                    if(varID2sample.containsKey(varID)){
                        varID2sample.get(varID).add(sample_head);
                    }else {
                        HashSet<String> newMap = new HashSet<>();
                        newMap.add(sample_head);
                        varID2sample.put(varID,newMap);
                    }

                    mergeVarMap.put(varID, varMap.get(varID));

                }
            }
        }


        String out_annofile = this.outdir + "/"+ prefix + this.output_var_table_suffix;
        String out_varfa = this.outdir + "/"+ prefix + this.output_db_table_suffix;

        BufferedWriter annoWriter = new BufferedWriter(new FileWriter(new File(out_annofile)));
        BufferedWriter pWriter = new BufferedWriter(new FileWriter(new File(out_varfa)));

        ArrayList<String> outSampleNames = new ArrayList<>();
        for(int i=0;i<sampleNames.size();i++){
            outSampleNames.add("sample:"+sampleNames.get(i));

        }
        annoWriter.write(JVariant.lineHead+"\t"+StringUtils.join(outSampleNames,"\t")+"\n");


        int k=0;
        for(String varID: mergeVarMap.keySet()){
            ArrayList<String> out = new ArrayList<>();
            for(String sample_name : sampleNames){
                ArrayList<String> ex = new ArrayList<>();
                for(int i=1;i<h.size();i++) {
                    String sample_head = sample_name + " " + h.get(i);
                    if(varID2sample.get(varID).contains(sample_head)){
                        ex.add("1");
                    }else{
                        ex.add("0");
                    }
                }

                out.add(StringUtils.join(ex,","));
            }

            String newVarID = "VAR|"+k;
            String d[] = mergeVarMap.get(varID).line.split("\t");
            d[0] = newVarID;
            annoWriter.write(StringUtils.join(d,"\t")+"\t"+StringUtils.join(out,"\t")+"\n");

            pWriter.write(">"+newVarID+" "+mergeVarMap.get(varID).protein_header+"\n"+mergeVarMap.get(varID).protein_seq+"\n");


            k++;
        }

        if(addRef){
            addReferenceProteins(gMap, pWriter, logWriter);
        }

        annoWriter.close();
        pWriter.close();

        String stat_file = this.outdir + "/"+ prefix + this.output_stat_table_suffix;
        printSummary(out_annofile,stat_file);


    }


    public HashMap<String, JVariant> readData(String annofile, String varfa, BufferedWriter logWriter) throws IOException {


        HashMap<String,String> proMap = new HashMap<>();
        HashMap<String,String> proInfoMap = new HashMap<>();

        File dbFile = new File(varfa);
        FASTAFileReader dbReader = new FASTAFileReaderImpl(dbFile);
        FASTAElementIterator it = dbReader.getIterator();
        // int num = 0;
        while (it.hasNext()) {
            FASTAElement el = it.next();
            el.setLineLength(1);
            String hLine = el.getHeader().trim();
            String headLine[] = hLine.split("\\s+");
            String pID = headLine[0];
            proMap.put(pID,el.getSequence());
            proInfoMap.put(pID,hLine.replaceAll("^[^ ]*",""));

        }
        dbReader.close();

        HashMap<String,JVariant> gMap = new HashMap<>();

        BufferedReader aReader = new BufferedReader(new FileReader(new File(annofile)));
        String headLine = aReader.readLine().trim();
        String h[] = headLine.split("\t");
        HashMap<String,Integer> hMap = new HashMap<>();
        for(int i=0;i<h.length;i++){
            hMap.put(h[i],i);
        }
        JVariant.lineHead = headLine;
        String line;
        while((line = aReader.readLine())!=null){

            line = line.trim();
            String d[] = line.split("\t");

            String varID = d[hMap.get("Chr")]+"_"+d[hMap.get("Start")]+"_"+d[hMap.get("End")]+"_"+d[hMap.get("Ref")]+"_"+d[hMap.get("Alt")]+"_"+d[hMap.get("mRNA")];
            JVariant jVariant = new JVariant();
            jVariant.ID = varID;
            jVariant.line = line;
            if(proMap.containsKey(d[hMap.get("Variant_ID")])) {
                jVariant.protein_seq = proMap.get(d[hMap.get("Variant_ID")]);
                jVariant.protein_header = proInfoMap.get(d[hMap.get("Variant_ID")]);
            }else{
                logWriter.write("Don't find protein sequence for variant: "+ d[hMap.get("Variant_ID")] + "("+line+") in file "+annofile+", "+varfa+"\n");
                continue;
            }

            gMap.put(varID,jVariant);

        }
        aReader.close();

        return(gMap);

    }


    private ArrayList<String> run(String variantAnno, String outfilePrefix, boolean addRef, boolean verbose, boolean addAnnoInfo, BufferedWriter logWriter) throws IOException {

        ////////////////////////////////////////////////////////////////////////////////////////////////////
        // output variant information
        String variant_table_file = this.outdir + "/" + outfilePrefix + this.output_var_table_suffix;
        String outdb = this.outdir + "/" + outfilePrefix + this.output_db_table_suffix;
        //String variant_table_file = this.outdb;
        //variant_table_file = variant_table_file.replaceAll(".fa$","_variant.txt");
        //if(variant_table_file.equalsIgnoreCase(this.outdb)){
        //    variant_table_file = variant_table_file.replaceAll(".fasta$","_variant.txt");
        //}

        //if(variant_table_file.equalsIgnoreCase(this.outdb)){
        //    variant_table_file = variant_table_file+"_variant.txt";
        //}

        BufferedWriter vWriter = new BufferedWriter(new FileWriter(new File(variant_table_file)));
        vWriter.write("Variant_ID\tChr\tStart\tEnd\tRef\tAlt\tVariant_Type\tVariant_Function\tGene\tmRNA\tProtein\tmRNA_Change\tProtein_Change\tAA_Ref\tAA_Pos\tAA_Var\n");

        ////////////////////////////////////////////////////////////////////////////////////////////////////
        // read variant annotation result by annovar
        BufferedReader vReader = new BufferedReader(new FileReader(new File(variantAnno)));
        // there is no header information
        String line ;
        Pattern p  = Pattern.compile("^[\\w\\-\\.]+:([\\w\\.]+):wholegene");
        Pattern p2 = Pattern.compile("^[\\w\\-\\.]+?:([\\w\\.]+?):exon\\d+:c.([\\w\\->]+)(:p.([\\w\\*]+))?");

        Pattern pcc1 = Pattern.compile("^([ACGTacgt])(\\d+)([ACGTacgt])$");
        Pattern pcc2 = Pattern.compile("^(\\d+)([ACGTacgt])>([ACGTacgt])$");
        // block substitution
        Pattern pcc3 = Pattern.compile("^(\\d+)_(\\d+)delins(\\w+)");
        // block substitution for a single nucleotide such as c.3delinsGT
        Pattern pcc4 = Pattern.compile("^(\\d+)delins(\\w+)");
        // single base deletion
        Pattern pcc5 = Pattern.compile("^(\\d+)del(\\w+)");
        // multi-base deletion
        Pattern pcc6 = Pattern.compile("^(\\d+)_(\\d+)del(\\w*)");
        // if end is equal to start, this is an insertion
        Pattern pcc7 = Pattern.compile("^(\\d+)_(\\d+)ins(\\w+)");
        // insertion
        Pattern pcc8 = Pattern.compile("^(\\d+)dup(\\w+)");
        // non-frameshift substitution (example: c.1825_1826TT, now the ref is not included in the string)
        Pattern pcc9 = Pattern.compile("^(\\d+)_(\\d+)(\\w+)");

        //
        HashMap<String,ArrayList<JVariant>> vMap = new HashMap<>();

        BufferedWriter pWriter = new BufferedWriter(new FileWriter(new File(outdb)));

        String headLine = vReader.readLine().trim();
        String hv[] = headLine.split("\t");
        HashMap<String,Integer> vheadMap = new HashMap<>();
        for(int i=0;i<hv.length;i++){
            vheadMap.put(hv[i],i);
        }

        HashMap<String,String> colNameMap = new HashMap<>();
        if(vheadMap.containsKey("AAChange.refGene")){

            colNameMap.put("Func","Func.refGene");
            colNameMap.put("Gene","Gene.refGene");
            colNameMap.put("ExonicFunc","ExonicFunc.refGene");
            colNameMap.put("AAChange","AAChange.refGene");

        }else if(vheadMap.containsKey("AAChange.ensGene")){

            colNameMap.put("Func","Func.ensGene");
            colNameMap.put("Gene","Gene.ensGene");
            colNameMap.put("ExonicFunc","ExonicFunc.ensGene");
            colNameMap.put("AAChange","AAChange.ensGene");

        }else if(vheadMap.containsKey("AAChange.knownGene")){

            colNameMap.put("Func","Func.knownGene");
            colNameMap.put("Gene","Gene.knownGene");
            colNameMap.put("ExonicFunc","ExonicFunc.knownGene");
            colNameMap.put("AAChange","AAChange.knownGene");

        }else{
            // for example: Gene.wgEncodeGencodeBasicV28
            // AAChange.wgEncodeGencodeBasicV28
            for(String i: vheadMap.keySet()){
                if(i.startsWith("AAChange")){
                    String suffix_name = i.replaceFirst("AAChange.","");
                    System.out.println("Annotation data:" + suffix_name);
                    colNameMap.put("Func","Func."+suffix_name);
                    colNameMap.put("Gene","Gene."+suffix_name);
                    colNameMap.put("ExonicFunc","ExonicFunc."+suffix_name);
                    colNameMap.put("AAChange","AAChange."+suffix_name);

                }
            }
        }

        // Func.refGene
        HashSet<String> funcMap = new HashSet<>();
        funcMap.add("exonic");
        funcMap.add("exonic;splicing");
        funcMap.add("splicing;exonic");


        int n_variant_id = 0;
        while((line = vReader.readLine())!=null){
            line = line.trim();

            String d[] = line.split("\t");

            String func_gene = d[vheadMap.get(colNameMap.get("Func"))];
            // exonic;splicing
            // exonic
            // ncRNA_exonic
            if(!funcMap.contains(func_gene)){
                continue;
            }

            String AAChange_Gene = d[vheadMap.get(colNameMap.get("AAChange"))];
            if(AAChange_Gene.equalsIgnoreCase("unknown")){
                continue;
            }else if(AAChange_Gene.equalsIgnoreCase(".")){
                continue;
            }

            // AAChange.refGene
            // 	$field[2] =~ m/^[\w\-\.]+:([\w\.]+):wholegene/ and next;
            Matcher m = p.matcher(AAChange_Gene);
            if(m.find()){
                continue;
            }

            // synonymous SNV removed
            String ExonicFunc_Gene = d[vheadMap.get(colNameMap.get("ExonicFunc"))];
            if(ExonicFunc_Gene.equalsIgnoreCase("synonymous SNV")){
                continue;
            }

            Matcher m2 = p2.matcher(AAChange_Gene);
            if(!m2.find()){
                System.err.println("Error: invalid record found in exonic_variant_function file (exonic format error): <"+line+">");
                logWriter.write("Error: invalid record found in exonic_variant_function file (exonic format error): <"+line+">"+"\n");
                System.exit(1);
            }

            String chr = d[vheadMap.get("Chr")];

            int genome_start = Integer.valueOf(d[vheadMap.get("Start")]);
            int genome_end = Integer.valueOf(d[vheadMap.get("End")]);
            String variant_type = ExonicFunc_Gene;

            String genome_ref = d[vheadMap.get("Ref")];
            String genome_var = d[vheadMap.get("Alt")];



            // PERM1:NM_001291366:exon2:c.T1183C:p.L395L,PERM1:NM_001291367:exon3:c.T901C:p.L301L
            String md[] = AAChange_Gene.split(",");
            for(int i=0;i<md.length;i++){
                Matcher m3 = p2.matcher(md[i]);
                if(m3.find()){
                    String dd[] = md[i].split(":");
                    String geneID = dd[0];
                    String transcript = dd[1];
                    String cchange = dd[3].replaceAll("^c.","");

                    // It's possible that there is no p.L301L part
                    // such as: CDK11A:NM_001313896:exon18:c.1973_1980CAGTCAAG
                    String pchange = "";
                    if(dd.length>=5) {
                        pchange = dd[4].replaceAll("^p.","");
                    }

                    Matcher mc1 = pcc1.matcher(cchange);
                    Matcher mc2 = pcc2.matcher(cchange);
                    Matcher mc3 = pcc3.matcher(cchange);
                    Matcher mc4 = pcc4.matcher(cchange);
                    Matcher mc5 = pcc5.matcher(cchange);
                    Matcher mc6 = pcc6.matcher(cchange);
                    Matcher mc7 = pcc7.matcher(cchange);
                    Matcher mc8 = pcc8.matcher(cchange);
                    Matcher mc9 = pcc9.matcher(cchange);

                    int start = -100;
                    int end = -100;
                    String ref = "";
                    String obs = "";

                    if(mc1.find()){
                        start = Integer.valueOf(mc1.group(2));
                        end   = start;
                        ref   = mc1.group(1);
                        obs   = mc1.group(3);
                    }else if(mc2.find()){
                        start = Integer.valueOf(mc2.group(1));
                        end   = start;
                        ref   = mc2.group(2);
                        obs   = mc2.group(3);
                    }else if(mc3.find()){
                        start = Integer.valueOf(mc3.group(1));
                        end   = Integer.valueOf(mc3.group(2));
                        obs   = mc3.group(3);
                    }else if(mc4.find()){
                        start = Integer.valueOf(mc4.group(1));
                        end   = start;
                        ref   = "REF";
                        obs   = mc4.group(2);
                    }else if(mc5.find()){
                        start = Integer.valueOf(mc5.group(1));
                        end   = start;
                        ref   = mc5.group(2);
                        obs   = "";
                    }else if(mc6.find()){
                        start = Integer.valueOf(mc6.group(1));
                        end   = Integer.valueOf(mc6.group(2));
                        ref   = mc6.group(3);
                        obs   = "";
                    }else if(mc7.find()){
                        start = Integer.valueOf(mc7.group(1));
                        end   = start;
                        ref   = "";
                        obs   = mc7.group(3);
                    }else if(mc8.find()){
                        start = Integer.valueOf(mc8.group(1));
                        end   = start;
                        ref   = "";
                        obs   = mc8.group(2);
                    }else if(mc9.find()){
                        start = Integer.valueOf(mc9.group(1));
                        end   = Integer.valueOf(mc9.group(2));
                        ref   = "REF";
                        obs   = mc9.group(3);
                    }else{
                        logWriter.write("Warning: invalid coding change format: <"+cchange+">"+ " within <"+ md[i]+"\n");
                        continue;
                    }

                    // variant

                    JVariant jVariant = new JVariant();
                    n_variant_id++;
                    jVariant.ID = String.valueOf(n_variant_id);// d[0];
                    jVariant.transcriptID = transcript;
                    jVariant.start = start;
                    jVariant.end = end;
                    jVariant.ref = ref;
                    jVariant.obs = obs;
                    jVariant.cchange = cchange;
                    jVariant.pchange = pchange;
                    jVariant.chr = chr;
                    jVariant.genome_start = genome_start;
                    jVariant.genome_end = genome_end;
                    jVariant.variant_type = variant_type;
                    jVariant.geneID = geneID;

                    jVariant.genome_ref = genome_ref;
                    jVariant.genome_var = genome_var;

                    if(vMap.containsKey(jVariant.ID)){
                        vMap.get(jVariant.ID).add(jVariant);
                    }else{
                        ArrayList<JVariant> js = new ArrayList<>();
                        js.add(jVariant);
                        vMap.put(jVariant.ID,js);
                    }



                }
            }


        }

        vReader.close();





        // process variants
        HashMap<String,Integer>  flagged_transcript = new HashMap<>();

        System.out.println("Transcripts with variant: "+vMap.size());
        logWriter.write("Transcripts with variant: "+vMap.size()+"\n");

        int nvar = 0;

        int cpu = Runtime.getRuntime().availableProcessors();

        ExecutorService fixedThreadPool = Executors.newFixedThreadPool(cpu);

        TranslateVariantWorker.gMap = gMap;


        for(String vID: vMap.keySet()){
            for(JVariant jv : vMap.get(vID)){
                nvar++;
                fixedThreadPool.execute(new TranslateVariantWorker(jv,nvar,verbose,logWriter));
            }
        }

        fixedThreadPool.shutdown();

        try {
            fixedThreadPool.awaitTermination(Long.MAX_VALUE, TimeUnit.HOURS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        for(String vID: vMap.keySet()){
            for(JVariant jv : vMap.get(vID)){
                if(jv.valid) {
                    vWriter.write(jv.line + "\n");

                    // it's possible to have a variant in which the mutant protein sequence is zero.
                    // In this case, we cannot output this item because it will cause problem when reading the fasta
                    if(jv.protein_seq != null && jv.protein_seq.length()>=1) {
                        pWriter.write(">" + jv.ID + " " + jv.protein_header + "\n");
                        pWriter.write(jv.protein_seq + "\n");
                    }else{
                        System.err.println("Protein length is zero: "+jv.ID);
                        logWriter.write("Protein length is zero: "+jv.ID+"\n");
                    }
                }
            }
        }




        // Add the reference protein sequences into the output database
        if(addRef){
            addReferenceProteins(gMap,pWriter,logWriter);
        }

        pWriter.close();
        vWriter.close();

        String stat_file = this.outdir + "/" + outfilePrefix + this.output_stat_table_suffix;
        printSummary(variant_table_file,stat_file);

        ArrayList<String> resfiles = new ArrayList<>();
        resfiles.add(variant_table_file);
        resfiles.add(outdb);
        return(resfiles);

    }

    public void readAnnotationData(String db, String geneAnno, BufferedWriter logWriter) throws IOException {
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        // read mRNA DNA sequences: hg19_refGeneMrna.fa
        File dbFile = new File(db);
        FASTAFileReader dbReader = new FASTAFileReaderImpl(dbFile);
        FASTAElementIterator it = dbReader.getIterator();
        int num = 0;
        //HashMap<String,JTranscript> gMap = new HashMap<>();

        // (leftmost exon at chr1:222001007)
        Pattern gInforPattern = Pattern.compile("leftmost exon at (.+?:\\d+)");


        HashSet<String> no_correct_orf_IDs = new HashSet<>();

        while (it.hasNext()) {
            FASTAElement el = it.next();
            el.setLineLength(1);
            String hLine = el.getHeader().trim();
            String headLine[] = hLine.split("\\s+");
            String pID = headLine[0];
            // some transcripts (such as NM_001075) occur multiple times in a file, sometimes with bad ORF annotation
            if(hLine.contains("does not have correct ORF annotation")){
                no_correct_orf_IDs.add(pID);
                continue;
            }


            Matcher gInforMatcher = gInforPattern.matcher(hLine);
            String infor = "";
            if(gInforMatcher.find()){
                infor = gInforMatcher.group(1);
            }

            String pSeq = el.getSequence().toUpperCase();

            if(gMap.containsKey(pID)){
                if(gMap.get(pID).mRNA_seq.length() < pSeq.length()){
                    gMap.get(pID).mRNA_seq = pSeq;
                    gMap.get(pID).infor = infor;
                }
            }else{
                JTranscript jTranscript = new JTranscript();
                jTranscript.transcriptID = pID;
                jTranscript.mRNA_seq = pSeq;
                jTranscript.infor = infor;
                gMap.put(pID,jTranscript);
            }
            num++;

        }
        System.out.println("Read mRNA sequences: "+num);
        logWriter.write("Read mRNA sequences: "+num+"\n");
        it.close();
        dbReader.close();

        String mRNA_used_fasta_file = this.outdir + "/mRNA_seq.fasta";
        String mRNA_used_anno_file = this.outdir + "/mRNA_anno.txt";
        System.out.println("Used mRNA sequences: "+mRNA_used_fasta_file);
        BufferedWriter mRNA_Seq_W = new BufferedWriter(new FileWriter(new File(mRNA_used_fasta_file)));
        for(String pID: gMap.keySet()){
            mRNA_Seq_W.write(">"+pID+" "+gMap.get(pID).infor+"\n"+gMap.get(pID).mRNA_seq+"\n");
        }
        mRNA_Seq_W.close();


        ////////////////////////////////////////////////////////////////////////////////////////////////////
        // read gene information: hg19_refGene.txt
        // no head
        // 585     NR_024540       chr1    -       14361   29370   29370   29370   11      14361,14969,15795,16606,16857,17232,17605,17914,18267,24737,29320,      14829,15038,15947,16765,17055,17368,17742,18061,18366,24891,29370,      0       WASH7P  unk     unk     -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        BufferedReader gReader = new BufferedReader(new FileReader(new File(geneAnno)));

        // refGene and ensGene has bin as the first column
        Pattern pGene = Pattern.compile("^\\d+$");

        // delete any transcripts that occurs multiple times in human genome
        HashMap<String, Integer> transcriptCopy = new HashMap<>();
        String line;
        while((line = gReader.readLine())!=null){
            line = line.trim();
            String d[] = line.split("\t");
            if(d.length<11){
                System.err.println("Error: invalid record found in gene file (>=11 fields expected): <"+line+">");
                logWriter.write("Error: invalid record found in gene file (>=11 fields expected): <"+line+">"+"\n");
                System.exit(1);
            }



            String name;
            String chr;
            String strand;
            int txstart;
            int txend;
            int cdsstart;
            int cdsend;
            String exonstart;
            String exonend;

            //
            Matcher pMatch = pGene.matcher(d[0]);
            if(pMatch.find()) { // refGene and ensGene has bin as the first column
                // remove the first element
                name = d[1];
                chr = d[2];
                strand = d[3];
                txstart = Integer.valueOf(d[4]);
                txend = Integer.valueOf(d[5]);
                cdsstart = Integer.valueOf(d[6]);
                cdsend = Integer.valueOf(d[7]);
                exonstart = d[9];
                exonend = d[10];
            }else{
                name = d[0];
                chr = d[1];
                strand = d[2];
                txstart = Integer.valueOf(d[3]);
                txend = Integer.valueOf(d[4]);
                cdsstart = Integer.valueOf(d[5]);
                cdsend = Integer.valueOf(d[6]);
                exonstart = d[8];
                exonend = d[9];
            }

            if(transcriptCopy.containsKey(name)){
                transcriptCopy.put(name,transcriptCopy.get(name)+1);
            }else{
                transcriptCopy.put(name,1);
            }

            if(!gMap.containsKey(name)){
                if(no_correct_orf_IDs.contains(name)) {
                    logWriter.write("Warning: remove (" + name +":does not have correct ORF annotation" + ") " + line + "\n");
                }else{
                    logWriter.write("Warning: remove (" + name + ") " + line + "\n");
                }
                continue;
            }

            if(!gMap.get(name).infor.isEmpty()){
                String infor = chr+":"+exonstart.split(",")[0];
                if(!gMap.get(name).infor.equalsIgnoreCase(infor)){
                    logWriter.write("Warning: remove (used "+name+" "+gMap.get(name).infor+") "+line+"\n");
                    continue;
                }
            }


            // we must save the original cds position information. we can know which gene is coding gene based on
            // these information.
            gMap.get(name).cds_start = cdsstart;
            gMap.get(name).cds_end = cdsend;


            // next we need to make sure that there is no intron between transcription start and translation start (this is rare but it happens when cdsstart is not in the first exon)
            exonstart = exonstart.replaceAll(",$","");
            exonend = exonend.replaceAll(",$","");
            String sexonstart[] = exonstart.split(",");
            String sexonend[] = exonend.split(",");

            // System.out.println("exonstart:"+exonstart+"===exonend:"+exonend);

            int dexonstart[] = new int[sexonstart.length];
            int dexonend[] = new int[sexonend.length];
            for(int i=0;i<sexonstart.length;i++){
                dexonstart[i] = Integer.valueOf(sexonstart[i]) + 1;
            }

            for(int i=0;i<sexonend.length;i++){
                dexonend[i] = Integer.valueOf(sexonend[i]);
            }

            txstart = txstart + 1;
            cdsstart = cdsstart + 1;

            int mrnastart = 0;
            int mrnaend = 0;

            if(strand.equalsIgnoreCase("+")){
                int intron = 0;
                for(int i=0;i<dexonstart.length;i++){
                    if(i!=0){
                        intron = intron + dexonstart[i] - dexonend[i-1] - 1;
                    }

                    if(cdsstart >= dexonstart[i] && cdsstart <= dexonend[i]){
                        mrnastart = cdsstart - txstart + 1 - intron;
                    }


                    if(cdsend >= dexonstart[i] && cdsend <= dexonend[i]){
                        mrnaend = cdsend - txstart + 1 - intron;
                    }

                    //System.out.println("i="+i+","+mrnastart+","+mrnaend);

                }

            }else if(strand.equalsIgnoreCase("-")){
                int intron = 0;
                for(int i=dexonstart.length-1;i>=0;i--){
                    if(i < (dexonstart.length-1)){
                        intron = intron + dexonstart[i+1] - dexonend[i] - 1;
                    }
                    if(cdsend >= dexonstart[i] && cdsend <= dexonend[i]){
                        mrnastart = txend - cdsend + 1 - intron;
                    }
                    if(cdsstart >= dexonstart[i] && cdsstart <= dexonend[i]){
                        mrnaend = txend - cdsstart + 1 - intron;
                    }
                }
            }

            gMap.get(name).start = mrnastart;
            gMap.get(name).end = mrnaend;
            gMap.get(name).transcriptID = name;


            // System.out.println(name+"\t"+mrnastart+"\t"+mrnaend);
            // gMap.put(name,jTranscript);

        }

        gReader.close();

        // note:
        // delete any transcripts that occurs multiple times in human genome
        int n_remove_transcript = 0;
        for(String tname : transcriptCopy.keySet()){
            if(transcriptCopy.get(tname) >=2){
                gMap.remove(tname);
                n_remove_transcript++;
            }
        }
        System.out.println("Remove transcripts: "+n_remove_transcript);
        logWriter.write("Remove transcripts: "+n_remove_transcript+"\n");

    }

    public void addReferenceProteins(HashMap<String,JTranscript> gMap, BufferedWriter pWriter, BufferedWriter logWriter) throws IOException {

        Pattern pStop = Pattern.compile("\\*.+");
        int nrefProteins = 0;

        for(String transcript: gMap.keySet()){
            if(gMap.containsKey(transcript)){
                // String dna_seq = mrnaMap.get(transcript).substring(gMap.get(transcript).start-1, gMap.get(transcript).end);

                // remove non-coding genes;
                if(gMap.get(transcript).cds_start == gMap.get(transcript).cds_end){
                    continue;
                }
                // System.out.println(transcript+"\t"+gMap.get(transcript).start+"\t"+gMap.get(transcript).end);
				String dna_seq = "";
				try {
                	dna_seq = gMap.get(transcript).mRNA_seq.substring(gMap.get(transcript).start-1, gMap.get(transcript).end);
				} catch (StringIndexOutOfBoundsException e) {
					e.printStackTrace();
                    logWriter.write(transcript+"\t"+gMap.get(transcript).mRNA_seq.length()+"\t"+gMap.get(transcript).start+"\t"+gMap.get(transcript).end+"\n");
					System.exit(1);
				}
                String protein = translateDNA(dna_seq,logWriter);
                Matcher pStopMatcher = pStop.matcher(protein);
                if(pStopMatcher.find()){
                    logWriter.write("Warning:"+transcript+" "+gMap.get(transcript).cds_start+" "+ gMap.get(transcript).cds_end+"\t"+protein+"\n");
                }

                protein = protein.replaceAll("\\*.*$","");
                try {
                    pWriter.write(">"+transcript+"\n"+protein+"\n");
                } catch (IOException e) {
                    e.printStackTrace();
                }
                nrefProteins++;


            }

        }
        System.out.println("Reference proteins: "+nrefProteins);
        logWriter.write("Reference proteins: "+nrefProteins+"\n");

    }

    public void addReferenceProteins(HashMap<String,JTranscript> gMap, String file, BufferedWriter logWriter) throws IOException {
        BufferedWriter pWriter = new BufferedWriter(new FileWriter(new File(file)));
        addReferenceProteins(gMap,pWriter,logWriter);
        pWriter.close();
    }


    public static String translateDNA(String dna,BufferedWriter logWriter) throws IOException {
        String protein = "";
        dna = dna.toUpperCase();
        Pattern p3 = Pattern.compile("(...)");
        Matcher matcher = p3.matcher(dna);
        while(matcher.find()){
            String t3 = matcher.group(1);
            if(codons.containsKey(t3)){
                protein = protein+codons.get(t3);
            }else{
                System.err.println("Error: invalid triplets found in DNA sequence to be translated: <"+t3+"> in <"+dna+">");
                logWriter.write("Error: invalid triplets found in DNA sequence to be translated: <"+t3+"> in <"+dna+">"+"\n");
                System.exit(1);
            }
        }

        return(protein);
    }

    public HashMap<String,String> getCodon(){

        HashMap<String,String> codonMap = new HashMap<>();
        codonMap.put("TTT","F");
        codonMap.put("TTC","F");
        codonMap.put("TCT","S");
        codonMap.put("TCC","S");
        codonMap.put("TAT","Y");
        codonMap.put("TAC","Y");
        codonMap.put("TGT","C");
        codonMap.put("TGC","C");
        codonMap.put("TTA","L");
        codonMap.put("TCA","S");
        codonMap.put("TAA","*");
        codonMap.put("TGA","*");

        codonMap.put("TTG","L");
        codonMap.put("TCG","S");
        codonMap.put("TAG","*");
        codonMap.put("TGG","W");
        codonMap.put("CTT","L");
        codonMap.put("CTC","L");
        codonMap.put("CCT","P");
        codonMap.put("CCC","P");
        codonMap.put("CAT","H");
        codonMap.put("CAC","H");
        codonMap.put("CGT","R");
        codonMap.put("CGC","R");
        codonMap.put("CTA","L");
        codonMap.put("CTG","L");
        codonMap.put("CCA","P");
        codonMap.put("CCG","P");
        codonMap.put("CAA","Q");
        codonMap.put("CAG","Q");
        codonMap.put("CGA","R");
        codonMap.put("CGG","R");
        codonMap.put("ATT","I");
        codonMap.put("ATC","I");
        codonMap.put("ACT","T");
        codonMap.put("ACC","T");
        codonMap.put("AAT","N");
        codonMap.put("AAC","N");
        codonMap.put("AGT","S");
        codonMap.put("AGC","S");
        codonMap.put("ATA","I");
        codonMap.put("ACA","T");
        codonMap.put("AAA","K");
        codonMap.put("AGA","R");
        codonMap.put("ATG","M");
        codonMap.put("ACG","T");
        codonMap.put("AAG","K");
        codonMap.put("AGG","R");
        codonMap.put("GTT","V");
        codonMap.put("GTC","V");
        codonMap.put("GCT","A");
        codonMap.put("GCC","A");
        codonMap.put("GAT","D");
        codonMap.put("GAC","D");
        codonMap.put("GGT","G");
        codonMap.put("GGC","G");
        codonMap.put("GTA","V");
        codonMap.put("GTG","V");
        codonMap.put("GCA","A");
        codonMap.put("GCG","A");
        codonMap.put("GAA","E");
        codonMap.put("GAG","E");
        codonMap.put("GGA","G");
        codonMap.put("GGG","G");
        return(codonMap);
    }


    /**
     *
     * @param anno_table_file The input is *-varInfo.txt
     */
    public void printSummary(String anno_table_file, String stat_file) throws IOException {


        //
        HashMap<String,HashSet<String>> typeMap = new HashMap<>();

        HashMap<String,Integer> typeMap_mRNA = new HashMap<>();

        HashMap<String,HashSet<String>> funMap = new HashMap<>();

        HashMap<String,Integer> funMap_mRNA = new HashMap<>();

        BufferedReader aReader = new BufferedReader(new FileReader(new File(anno_table_file)));

        String headLine = aReader.readLine().trim();
        String h[] = headLine.split("\t");
        HashMap<String,Integer> hMap = new HashMap<>();
        for(int i=0;i<h.length;i++){
            hMap.put(h[i],i);
        }

        String line;

        // They are used to save all the variant types and functions for output.
        HashSet<String> allTypes = new HashSet<>();
        HashSet<String> allFunctions = new HashSet<>();

        HashSet<String> all_variant_no_mRNA = new HashSet<>();

        int total_variants = 0;

        while((line = aReader.readLine())!=null){
            String d[] = line.trim().split("\t");
            String variant_type = d[hMap.get("Variant_Type")];
            String variant_fun = d[hMap.get("Variant_Function")];

            allTypes.add(variant_type);
            allFunctions.add(variant_fun);

            // don't consider mRNA
            String varID = d[hMap.get("Chr")]+"_"+d[hMap.get("Start")]+"_"+d[hMap.get("End")]+"_"+d[hMap.get("Ref")]+"_"+d[hMap.get("Alt")];

            all_variant_no_mRNA.add(varID);
            // variant type
            if(typeMap.containsKey(variant_type)){
                typeMap.get(variant_type).add(varID);
            }else{
                HashSet<String> vlist = new HashSet<>();
                vlist.add(varID);
                typeMap.put(variant_type,vlist);
            }

            // function type
            if(funMap.containsKey(variant_fun)){
                funMap.get(variant_fun).add(varID);
            }else{
                HashSet<String> flist = new HashSet<>();
                flist.add(varID);
                funMap.put(variant_fun,flist);
            }

            // consider mRNA
            // varID = d[hMap.get("Chr")]+"_"+d[hMap.get("Start")]+"_"+d[hMap.get("End")]+"_"+d[hMap.get("Ref")]+"_"+d[hMap.get("Alt")]+"_"+d[hMap.get("mRNA")];
            if(typeMap_mRNA.containsKey(variant_type)){
                typeMap_mRNA.put(variant_type,typeMap_mRNA.get(variant_type)+1);
            }else{
                typeMap_mRNA.put(variant_type,1);
            }
            if(funMap_mRNA.containsKey(variant_fun)){
                funMap_mRNA.put(variant_fun,funMap_mRNA.get(variant_fun)+1);
            }else{
                funMap_mRNA.put(variant_fun,1);
            }

            total_variants++;


        }

        aReader.close();

        BufferedWriter pWriter = new BufferedWriter(new FileWriter(new File(stat_file)));

        pWriter.write("#Variant_Type\n");
        for(String type: allTypes){
            double ratio1 = 1.0*typeMap.get(type).size()/all_variant_no_mRNA.size();
            double ratio2 = 1.0*typeMap_mRNA.get(type)/total_variants;
            pWriter.write(type+"\t"+typeMap.get(type).size()+"\t"+ratio1+"\t"+typeMap_mRNA.get(type)+"\t"+ratio2+"\n");
        }

        pWriter.write("#Variant_Function\n");
        for(String fun: allFunctions){
            double ratio1 = 1.0*funMap.get(fun).size()/all_variant_no_mRNA.size();
            double ratio2 = 1.0*funMap_mRNA.get(fun)/total_variants;
            pWriter.write(fun+"\t"+funMap.get(fun).size()+"\t"+ratio1+"\t"+funMap_mRNA.get(fun)+"\t"+ratio2+"\n");
        }


        pWriter.close();



    }



















}
