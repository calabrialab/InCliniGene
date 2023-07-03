import java.nio.file.*;

import java.io.IOException;
import java.util.*;

import java.io.FileWriter;
import java.io.BufferedWriter;

import java.io.BufferedReader;
import java.io.FileReader;

import java.util.Collections;



class Insertion {
    String clu_id = null;
    String cons_seq = null;

    //int num_aln = 1;
    //int max_aln_score = 1;
    //int seq_len = 1;
    int weight;

    String target_id = null;
    int centroid;
    int offset;
    int aln_score = 1;

    public Insertion(int id, String cs, int w, String tid, int c, int o) {
        clu_id = "subg_"+id;
        cons_seq = cs;
        weight = w;
        target_id = "chr"+tid;
        centroid = c;
        offset = o;
    }

    /*
 ISCluster ::= {
  clu_id "subg_1", //valore incrementale
  cons_seq "cippa", /valore arbitrario, VB nel seguito
  num_aln 1,  //VB
  max_aln_score -1, //VB 
  seq_len -1, //VB
  weight 61,
  clu_score { 1, 10, 0 }, //VB
  lab {
    {
      target_id "chr10",
      centroid { 100026030, 10, 0 },
      offset +1, //avrò solo +1 e -1
      aln_score -1, //VB
      seq_id {
        "subg_1" //Stesso nome di sopra, ma tanto questo non lo leggo
      }
    }
  },
  merged {
  },
  masked {
  }
}

 */
    public void printInsertion(BufferedWriter logfile) {
        try {
            logfile.write(" ISCluster ::= {\n");
            logfile.write("  clu_id \""+clu_id+"\",\n");
            logfile.write("  cons_seq \""+cons_seq+"\",\n");
            logfile.write("  num_aln 1,\n");
            logfile.write("  max_aln_score 1,\n");
            logfile.write("  seq_len 1,\n");
            logfile.write("  weight "+weight+",\n");
            logfile.write("  clu_score { 1, 10, 0 },\n");
            logfile.write("  lab {\n");
            logfile.write("    {\n");
            logfile.write("      target_id \""+target_id+"\",\n");
            logfile.write("      centroid { "+centroid+", 10, 0 },\n");
            logfile.write("      offset "+offset+", \n");
            logfile.write("      aln_score 1,\n");
            logfile.write("      seq_id {\n");
            logfile.write("        \""+clu_id+"\"\n");
            logfile.write("      }\n");
            logfile.write("   }\n");
            logfile.write("  },\n");
            logfile.write("  merged {\n");
            logfile.write("  },\n");
            logfile.write("  masked {\n");
            logfile.write("  },\n");
            logfile.write("}\n");

        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}



class Sample {
    int numsubg = 1;
    String name;
    Vector<Insertion> insertionlist = new Vector<Insertion>(1);

    public Sample(String n) {
        name = n;
    }

    public void addInsertion(String cs, int w, String tid, int c, int o) {
        insertionlist.addElement(new Insertion(numsubg, cs, w, tid, c, o));
        numsubg++;
    }

    public void printSample() {
        BufferedWriter outfile;

        try {
            outfile = new BufferedWriter(new FileWriter(name+".clu"));

            Enumeration en = insertionlist.elements();

            while(en.hasMoreElements()) {
                Insertion i = (Insertion)en.nextElement();
                i.printInsertion(outfile);
            }
            outfile.close();


        }  catch (IOException x) {
            System.err.println(x);
        }




    }
}
    



public class sparseMatrix2Jason {


    public static void main(String[] args) {
        Path file = Paths.get("./AssayValidation-CLONALEXPANSION-allPools_seqCount_matrix.no0.annotated.tsv");
        // Path file = null;
        // Path file = StdIn.readLine();

        boolean gs = false;
        String chr;	
        String  integration_locus = "0";
        int strand = 0;
        String GeneName ="NN", GeneStrand;

        Vector<Sample> vsample = new Vector<Sample>(1);

        try {
            BufferedReader br = new BufferedReader(new FileReader(file.toString()));
            String line, lineaux;
            int num_insertion_file = 0;
            int numbrackets = 0;
            boolean isfirst = true;

            int ncol = 0, samplenum = 0;

            String[] wordsArray;

            line = br.readLine();
 
            wordsArray = line.split("\t");
            for(String each : wordsArray){
                if(!"".equals(each)){
                    if("chr".equals(each)) continue;
                    if("integration_locus".equals(each)) continue;
                    if("strand".equals(each)) continue;
                    if("GeneName".equals(each)) { gs = true; continue; }
                    if("GeneStrand".equals(each)) continue;
 
                    vsample.addElement(new Sample(each));
                }
            }

            while ((line = br.readLine()) != null) {
                wordsArray = line.split("\t");
                ncol = 0;
                samplenum = 0;
                chr="";
                for(String each : wordsArray){
                    if(!"".equals(each)){
                        if (ncol == 0) { chr = each; ncol++; continue; }	
                        if (ncol == 1) { integration_locus = each; ncol++; continue;}
                        if (ncol == 2) { 
                            if(each.equals("+")) strand = 1;
                            else strand = -1; 
                            ncol++; 
                            continue;
                        }
                        if (ncol == 3 && gs) { GeneName = each; ncol++; continue;}
                        else if (ncol == 3 && !gs) ncol = 5;
                        if (ncol == 4 && gs) { GeneStrand = each; ncol++; continue;}

                        vsample.get(samplenum).addInsertion(GeneName,Math.round(Float.parseFloat(each)),chr,
                            Integer.parseInt(integration_locus),strand);
                    }
                    if (ncol == 5) samplenum++;
                }
            }

            br.close();
        }
        catch (IOException x) {
            System.err.println(x);
        }

        Enumeration en = vsample.elements();

        while(en.hasMoreElements()) {
            Sample t = (Sample)en.nextElement();
            t.printSample();
        }

        System.out.println("DONE");
        System.exit(0);

    }
}