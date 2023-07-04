import org.neo4j.driver.v1.*;
import org.neo4j.driver.v1.exceptions.NoSuchRecordException;

import java.nio.file.*;

import java.io.IOException;
import java.util.*;

import java.io.FileWriter;
import java.io.BufferedWriter;

import java.io.BufferedReader;
import java.io.FileReader;

import java.util.Collections;



class Target implements Comparable<Target>{
    String target_id = null;
    int centroid = 0;
    int offset = 0;
    int aln_score = 0;

    public Target(String tid, int cent, int off, int as) {
        target_id = tid;
        centroid = cent;
        offset = off;
        aln_score = as;
    }

    public void PrintTarget() {
        System.out.println(target_id+"\t"+centroid+"\t"+offset+"\t"+aln_score);
    }

    public int compareTo(Target otherT) {
        return this.aln_score - otherT.aln_score;
    }
}


/******************************************************/
/******************************************************/
/******************************************************/



class Insertion {
    String clu_id = null;
    String cons_seq = null;
    Vector<Target> targetlist = new Vector<Target>(1);
    Integer num_aln = null;
    Integer max_aln_score = null;
    Integer seq_len = null;
    Integer weight = null;
    Vector<String> masked = new Vector<String>(1);


    public void PrintInsertion(String id) {
        Enumeration en = targetlist.elements();

        System.out.println("*** "+id+"_"+clu_id+ " T: "+targetlist.size() + " L: "+masked.size()+" ***");
        while(en.hasMoreElements()) {
            Target t = (Target)en.nextElement();
            t.PrintTarget();

        }
        if(!masked.isEmpty()) {
            //https://www.javacodeexamples.com/remove-duplicate-elements-from-vector-in-java-example/3271

            LinkedHashSet<String> lhSet = new LinkedHashSet<String>(masked);
            masked.clear();
            masked.addAll(lhSet);


            en = masked.elements();
            while(en.hasMoreElements()) {
                String s = (String) en.nextElement();
                System.out.println("--- "+s);
            }
        }

    }




    public void AddInsertion(String id, Driver d, BufferedWriter logfile, Hashtable<String, Integer> hsample) {
        String cluid = id + "_" + clu_id;
        String prefixclu = id + "_";

        Enumeration en;
        String command1, command2, command3;
        StatementResult result;
        String subgres;

        String candidate;
        int degree1, degree2;
        Integer count;

        long startTime, stopTime;
        Float g4val = Float.valueOf(0);

        Hashtable<String, Integer> hmpos = new Hashtable<String, Integer>();    //To count how many genomic positions are shared
        Hashtable<String, Integer> hmlabel = new Hashtable<String, Integer>();  //To check the label sharing
        Hashtable<String, Integer> hmdegree = new Hashtable<String, Integer>(); //To get the degree of the considered subg
        Hashtable<String, Float> hmmain = new Hashtable<String, Float>();

        System.out.println("Starting " + id + "_" + clu_id + " T : " + num_aln + " ***");

        //Collegamento tra sample e subg
        command1 = "MATCH (ss:sample {UniqueID: \"" + id + "\" }) \n";
        command1 = command1.concat("MERGE (s:subg {clu_id: \"" + cluid + "\", " +
                //"cons_seq: \"" + cons_seq + "\", "+
                "num_aln: " + num_aln + ", " +
                "cons_seq: \"" + cons_seq + "\", " +
                "max_aln_score: " + max_aln_score + ", " +
                "seq_len: " + seq_len + ", " +
                "weight: " + weight + "})\n");
        command1 = command1.concat("MERGE (ss)-[:sample2subg]->(s)");


        try (Session session = d.session()) {
            session.run(command1);
        }

        command1 = "MATCH (s:subg {clu_id: \"" + cluid + "\"})\n";
        command3 = null;

        /***********************************  POS ***********************************/
        //(n:subg)-[l:insertion]-(m:pos)
        int i = 0;
        en = targetlist.elements();

        startTime = System.currentTimeMillis(); ////////////////////////////

        command2 = "MATCH (s:subg {clu_id: $clid}) " +
                "MERGE (n:pos {target_id: $tid, centroid: $ci, searchlabel: $sl}) " +
                "MERGE (s)-[:insertion {offset: $os, aln_score: $as}]-(n)";

        try (Transaction tx = d.session().beginTransaction()) {
            while (en.hasMoreElements()) {
                Target t = (Target) en.nextElement();


                //Costruisco la ricerca di altre subg di sample diversi che condividono l'inserzione con intorno [+2,-2]
                if (command3 != null) command3 = command3.concat(", \"" + t.target_id + "_" + t.centroid + "\", \"" +
                        t.target_id + "_" + Integer.toString(t.centroid - 2) + "\", \"" + t.target_id + "_" + Integer.toString(t.centroid - 1) + "\", \"" +
                        t.target_id + "_" + Integer.toString(t.centroid + 1) + "\", \"" + t.target_id + "_" + Integer.toString(t.centroid + 2)) + "\" ";
                else
                    command3 = "MATCH (b:pos)-[:insertion]-(a:subg) where NOT a.clu_id STARTS WITH \"" + prefixclu + "\" and  " +
                            "b.searchlabel IN [\"" + t.target_id + "_" + t.centroid + "\", \"" +
                            t.target_id + "_" + Integer.toString(t.centroid - 2) + "\", \"" + t.target_id + "_" + Integer.toString(t.centroid - 1) + "\", \"" +
                            t.target_id + "_" + Integer.toString(t.centroid + 1) + "\", \"" + t.target_id + "_" + Integer.toString(t.centroid + 2) + "\" ";

                //System.out.println(command3);
                i++;

                tx.run(command2, Values.parameters("clid", cluid, "tid", t.target_id, "ci", t.centroid, "sl", t.target_id + "_" + t.centroid, "os", t.offset, "as", t.aln_score));
            }
            tx.success();
        }


        stopTime = System.currentTimeMillis(); ////////////////////////////
        System.out.println("Position took " + (stopTime - startTime) + " milliseconds");

        command3 = command3.concat("]  return a.clu_id, a.num_aln ");
        //System.out.println(command3);


        startTime = System.currentTimeMillis(); ////////////////////////////

        try (Session session = d.session()) {
            result = session.run(command3);

            while (result.hasNext()) {
                Record record = result.next();
                subgres = record.get(0).asString();
                count = hmpos.get(subgres);
                hmpos.put(subgres, count == null ? 1 : ++count);

                //Inserisco il degree della subg per determinare il tipo di match
                /*
                try {  ///DA MODIFICARE!
                    count = Integer.parseInt(record.get(1).asString());
                }  catch (Exception e) {
                    count = record.get(1).asInt();
                }
                hmdegree.put(subgres, count);*/
                hmdegree.put(subgres, record.get(1).asInt());
            }
        }


        stopTime = System.currentTimeMillis(); ////////////////////////////
        System.out.println("Query took " + (stopTime - startTime) + " milliseconds");


        degree1 = i;
        if (degree1 != num_aln) {
            System.out.println("*** *** *** Problema sul numero di allineamenti "+degree1+ " e "+ num_aln+" *** *** ***");
            try {
                logfile.write("*** *** *** Problema sul numero di allineamenti "+degree1+ " e "+ num_aln+" *** *** ***\n");
            } catch (IOException e) {
                e.printStackTrace();
            }
        }



        /* Search for MAIN INSERTION for case 4 */
        if (degree1 > 1) { // The two topmost insertions differs for at least 30%

            startTime = System.currentTimeMillis(); ////////////////////////////
            Collections.sort(targetlist, Collections.reverseOrder());
            stopTime = System.currentTimeMillis(); ////////////////////////////

            try {
                logfile.write("SORT took " + (stopTime - startTime) + " milliseconds\n");
            } catch (IOException e) {
                e.printStackTrace();
            }

            float aa = (float) targetlist.get(0).aln_score;
            float bb = (float) targetlist.get(1).aln_score;

            if ( (bb/aa) < 0.7 ) {
                g4val = bb/aa;
                String mtg = targetlist.get(0).target_id;
                int mco = targetlist.get(0).centroid;
                int mof = targetlist.get(0).offset;

                startTime = System.currentTimeMillis(); ////////////////////////////

                String commanda = "MATCH (s:subg {clu_id: $clid}) " +
                        "MERGE (n:pos {target_id: $tid, centroid: $ci}) " +
                        "MERGE (s)-[:main_insertion {offset: $os, aln_score: $as, val: $v}]-(n)";

                //SEARCH for other MAIN INSERTIONS on the same genomic position
                String commandb = "MATCH (b:pos)-[l:main_insertion]-(a:subg) where NOT a.clu_id STARTS WITH \"" + prefixclu +
                        "\" and b.searchlabel IN [\"" + mtg + "_" + mco + "\", \"" +
                        mtg + "_" + Integer.toString(mco - 1) + "\", \"" + mtg + "_" + Integer.toString(mco +1) + "\", \"" +
                        mtg + "_" + Integer.toString(mco - 2) + "\", \"" + mtg + "_" + Integer.toString(mco +2) +
                        "\" ]  return a.clu_id, l.val ";

                try (Session session = d.session()) {
                    session.run(commanda, Values.parameters("clid", cluid, "tid", mtg, "ci", mco, "os", mof, "as", aa, "v", g4val));

                    result = session.run(commandb);
                    while (result.hasNext()) {
                        Record record = result.next();
                        hmmain.put(record.get(0).asString(), record.get(1).asFloat());

                        try {
                            logfile.write("1 Found "+record.get(0).asString() +"\n");
                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                    }
                }

                //SEARCH for single INSERTION on the same genomic position
                commandb = "MATCH (b:pos)-[l:insertion]-(a:subg) where NOT a.clu_id STARTS WITH \"" + prefixclu +
                        "\" and a.num_aln = 1 and b.searchlabel IN [\"" + mtg + "_" + mco + "\", \"" +
                        mtg + "_" + Integer.toString(mco - 1) + "\", \"" + mtg + "_" + Integer.toString(mco +1) + "\", \"" +
                        mtg + "_" + Integer.toString(mco - 2) + "\", \"" + mtg + "_" + Integer.toString(mco +2) +"\" ]  return a.clu_id ";

                try (Session session = d.session()) {
                    result = session.run(commandb);
                    while (result.hasNext()) {
                        Record record = result.next();
                        //if (record.get(1).asInt() == 1 )
                        hmmain.put(record.get(0).asString(), Float.valueOf(-1));

                        try {
                            logfile.write("2 Found "+record.get(0).asString() +"\n");
                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                    }
                }

                stopTime = System.currentTimeMillis(); ////////////////////////////


                try {
                    logfile.write("*** MAIN_INSERTION "+mtg+ "_"+ mco+" : "+g4val+"\n");
                    //logfile.write(commandb+"\n");
                    logfile.write("It tooks " + (stopTime - startTime) + " milliseconds\n");
                } catch (IOException e) {
                    e.printStackTrace();
                }

            }
        } else { //A subg with a single assignemnt search for MAIN INSERTIONS
            String mtg = targetlist.get(0).target_id; //Only 1 pos
            int mco = targetlist.get(0).centroid;

            String commandb = "MATCH (b:pos)-[l:main_insertion]-(a:subg) where NOT a.clu_id STARTS WITH \"" + prefixclu +
                    "\" and b.searchlabel IN [\"" + mtg + "_" + mco + "\", \"" +
                    mtg + "_" + Integer.toString(mco - 1) + "\", \"" + mtg + "_" + Integer.toString(mco +1) + "\", \"" +
                    mtg + "_" + Integer.toString(mco - 2) + "\", \"" + mtg + "_" + Integer.toString(mco +2) +"\" ]  return a.clu_id, l.val ";

            try (Session session = d.session()) {
                result = session.run(commandb);
                while (result.hasNext()) {
                    Record record = result.next();
                    hmmain.put(record.get(0).asString(), record.get(1).asFloat());

                    try {
                        logfile.write("3 Found "+record.get(0).asString() +"\n");
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
            }
        }


        /***********************************  LABEL ***********************************/
        //(n:subg)-[l:repeat]-(m:label)

        command3 = null;
        i = 0;

        //SE ho etichette le inserisco
        if(!masked.isEmpty()) {

            startTime  = System.currentTimeMillis(); ////////////////////////////

            //https://www.javacodeexamples.com/remove-duplicate-elements-from-vector-in-java-example/3271

            LinkedHashSet<String> lhSet = new LinkedHashSet<String>(masked);
            masked.clear();
            masked.addAll(lhSet);

            en = masked.elements();
            while(en.hasMoreElements()) {
                String s = (String) en.nextElement();

                //Collegamento tra subg e label
                command2 = command1.concat("MERGE (m"+i+":label {name: \"" +s+ "\"})\n");
                command2 = command2.concat("MERGE (s)-[:repeat]->(m"+i+")\n");
                i++;

                //Altre subg che condividono la label
                if(command3 != null) command3 = command3.concat(", \""+s+"\"");
                else command3 = "MATCH (a:subg)-[:repeat]-(b:label) where NOT a.clu_id STARTS WITH \"" +prefixclu+ "\" and b.name IN  [\"" +s+ "\"";

                ///System.out.println(command3);

                try (Session session = d.session()) {
                    session.run(command2);
                }
            }


            command3 = command3.concat("] return a.clu_id ");
            //System.out.println(command3);

            try (Session session = d.session()) {
                result = session.run(command3);

                while ( result.hasNext() )
                {
                    //Qui non mi interessa il degree
                    subgres = result.next().get(0).asString();
                    count = hmlabel.get(subgres);
                    hmlabel.put(subgres, count == null ? 1 : ++count);
                }
            }

            stopTime = System.currentTimeMillis(); ////////////////////////////
            System.out.println("Label took " + (stopTime - startTime) + " milliseconds");
        }



        ////(n:subg)-[l:gtris]-(m:subg) m.prefixclu != n.prefixclu

        /* REGOLE per il collegamento
        1) se ha un'inserzione unica e l'altra pure a posizione +-3 OK (A2-B2-C2) NO LABEL
        2) se le etichette di allineamento (target, centroid) sono condivise almeno per il 50%  (e non ci sono etichette di repeat)
        3) condividono almeno una etichetta della repeat (almeno una) ed ALMENO UN punto di inserzione (etichetta/annotazione genomica) (sempre distanza +-3) (A3-C3, A1-B1)
        4) se ho inserzione unica ed una main_insertion o due main_insertion
        */



        startTime  = System.currentTimeMillis(); ////////////////////////////

        Set<String> keys;
        Iterator<String> itr;


        //Regole 1 e 2
        if (!hmpos.isEmpty()) {
            keys = hmpos.keySet();
            itr = keys.iterator();
            while (itr.hasNext()) {

                candidate = itr.next();
                degree2 = hmdegree.get(candidate);

                if(degree1 == 1 && degree2 == 1) { //REGOLA 1
                    command1 =  "MATCH (a:subg {clu_id: \"" +cluid+ "\"})\n"+
                                "MATCH (b:subg {clu_id: \"" +candidate+ "\"})\n"+
                                "MERGE (a)-[:gtris {case: \"1\"}]-(b)";
                    //command1 = ("MERGE (a:subg {clu_id: \"" +cluid+ "\"})-[:gtris {case: \"1\"}]-(b:subg {clu_id: \"" +candidate+ "\"})\n");

                    //Conto quanto link tra i due sample delle due subg
                    String[] candidateid = candidate.split("_");
                    count = hsample.get(candidateid[0]);
                    hsample.put(candidateid[0], count == null ? 1 : ++count);

                    hmlabel.remove(candidate); //Spesso 1 implica 3, per evitare di avere due link gtris
                    hmmain.remove(candidate); //NON dovrebbe capitare
                    try (Session session = d.session()) {
                        session.run(command1);
                        System.out.println("*** GTRIS 1 between "+cluid+ " and "+ candidate);
                        try {
                            logfile.write("*** GTRIS 1 between "+cluid+ " and "+ candidate);
                            logfile.write(" and between "+id+" and "+candidateid[0]+"\n");
                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                    }
                } else { //REGOLA 2
                    count = hmpos.get(candidate);

                    //SE OR avrei un GTRIS 2 monodirezionale
                    //Se un nodo ha solo 2 inserzioni considero 1 come il 50%+1, altrimenti dovrei avere esattamente 2

                    boolean left, right;
                    if (degree1 == 2 ) left = (count >= 1); //dovrebbe essere sempre vera
                    else left = (count > (degree1/2));

                    if (degree2 == 2 ) right = (count >= 1); //idem
                    else right = (count > (degree2/2));

                    //if (count > (degree1/2) && count > (degree2/2) ) {
                    if( left && right ) {
                        command1 =  "MATCH (a:subg {clu_id: \"" +cluid+ "\"})\n"+
                                    "MATCH (b:subg {clu_id: \"" +candidate+ "\"})\n"+
                                    "MERGE (a)-[:gtris {case: \"2\"}]-(b)";

                        String[] candidateid = candidate.split("_");
                        count = hsample.get(candidateid[0]);
                        hsample.put(candidateid[0], count == null ? 1 : ++count);

                        hmlabel.remove(candidate); //A volte
                        hmmain.remove(candidate);

                        try (Session session = d.session()) {
                            session.run(command1);
                            System.out.println("*** GTRIS 2 between "+cluid+ " and "+ candidate);
                            try {
                                logfile.write("*** GTRIS 2 between "+cluid+ " and "+ candidate);
                                logfile.write(" and between "+id+" and "+candidateid[0]+"\n");
                            } catch (IOException e) {
                                e.printStackTrace();
                            }
                        }
                    }

                }
            }
        }

        //Regola 3: due subg condividono almeno 1 label ed almeno 1 insertion
        if (!hmlabel.isEmpty()) {
            keys = hmlabel.keySet();
            itr = keys.iterator();

            while (itr.hasNext()) {  //Scorro tutte le subg che hanno in comune almeno una label
                candidate = itr.next();
                if (hmpos.get(candidate) != null) { //E' presente almeno una inserzione in comune
                    command1 =  "MATCH (a:subg {clu_id: \"" +cluid+ "\"})\n"+
                            "MATCH (b:subg {clu_id: \"" +candidate+ "\"})\n"+
                            "MERGE (a)-[:gtris {case: \"3\"}]-(b)";

                    //command1 = ("MERGE (a:subg {clu_id: \"" +cluid+ "\"})-[:gtris {case: \"3\"}]-(b:subg {clu_id: \"" +candidate+ "\"})\n");


                    //Conto quanti link tra i due sample delle due subg
                    String[] candidateid = candidate.split("_");
                    count = hsample.get(candidateid[0]);
                    hsample.put(candidateid[0], count == null ? 1 : ++count);

                    hmmain.remove(candidate); //potrebbe capitare

                    try (Session session = d.session()) {
                        session.run(command1);
                        //session.run(command2);
                        System.out.println("*** GTRIS 3 between "+cluid+ " and "+ candidate);
                        try {
                            logfile.write("*** GTRIS 3 between "+cluid+ " and "+ candidate );
                            logfile.write(" and between "+id+" and "+candidateid[0]+"\n");
                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                    }

                }
            }
        }

        //Regola 4:  inserzione unica ed una main_insertion o due main_insertion
        if (!hmmain.isEmpty()) {
            keys = hmmain.keySet();
            itr = keys.iterator();
            Float maxv;

            while (itr.hasNext()) {  //Scorro tutte le subg che hanno in comune almeno una label
                candidate = itr.next();

                if (g4val > hmmain.get(candidate)) maxv = g4val;
                else maxv = hmmain.get(candidate);

                command1 =  "MATCH (a:subg {clu_id: \"" +cluid+ "\"})\n"+
                        "MATCH (b:subg {clu_id: \"" +candidate+ "\"})\n"+
                        "MERGE (a)-[:gtris {case: \"4\", val: "+maxv+"}]-(b)";

                String[] candidateid = candidate.split("_");
                count = hsample.get(candidateid[0]);
                hsample.put(candidateid[0], count == null ? 1 : ++count);


                try (Session session = d.session()) {
                    session.run(command1);
                    //session.run(command2);
                    System.out.println("*** GTRIS 4 between "+cluid+ " and "+ candidate);
                    try {
                        logfile.write("*** GTRIS 4 between "+cluid+ " and "+ candidate );
                        logfile.write(" and between "+id+" and "+candidateid[0]+"\n");
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
            }
        }


        stopTime = System.currentTimeMillis(); ////////////////////////////
        System.out.println("gtris took " + (stopTime - startTime) + " milliseconds");

        //System.out.println(command);
        System.out.println("Added "+id+"_"+clu_id+ " T : "+targetlist.size() + " L : "+masked.size()+" ***");
        try {
            logfile.write("Added "+id+"_"+clu_id+ " T : "+targetlist.size() + " L : "+masked.size()+" ***\n");
        } catch (IOException e) {
            e.printStackTrace();
        }

    }
}



/******************************************************/
/******************************************************/
/******************************************************/




public class CluParser {


    public static void main(String[] args) {
        //Path file = Paths.get("/Users/dago/Desktop/ivan/BED/unique_r1_ID00000000000000016511.filt.out.all.noprobl.clu");
        Path file = null;
        String ID ="Error";

        Hashtable<String, Integer> hsample = new Hashtable<String, Integer>();

        Driver driver;

        if (file == null)
            if (args.length == 0) {
                System.out.println("Indicare il file da leggere! \n Uso: java CluParser filename [ID]");
                System.exit(-1);
            } else file = Paths.get(args[0]);
        
        //LOCAL ----- CHANGE HERE WITH YOUR UID AND PASSWORD (!!!!)
        driver = GraphDatabase.driver("bolt://localhost:7687", AuthTokens.basic(/*UID*/, /*PWD*/));
        
        if(args.length > 1) ID = args[1];
        else if (file.getFileName().toString().contains("unique_r1_ID")) {
            System.out.println("Il nome del file continene un ID: "); //E me lo cerco del posto opportuno

            int i1 = file.getFileName().toString().indexOf("_ID");
            int i2 = file.getFileName().toString().indexOf(".");

            ID = file.getFileName().toString().substring(i1+1,i2);
        } else if (file.getFileName().toString().contains("unique_r1_LMv") || file.getFileName().toString().contains("unique_r1_LAM")) {
            int i1 = file.getFileName().toString().indexOf("_r1");
            int i2 = file.getFileName().toString().indexOf("_LC");
            int i3 = file.getFileName().toString().indexOf(".filt");

            //DA LAM-LTR_Block_L_11_LTR21_LC66_BR2_TR1
            //A  LAM-LTR_Block_L_11_LTR21.LC66_BR2_TR1

            String CompleteAmplificationID = file.getFileName().toString().substring(i1+4,i2);
            CompleteAmplificationID = CompleteAmplificationID.concat("."+file.getFileName().toString().substring(i2+1,i3));

            String command1 = "MATCH (ss:sample {CompleteAmplificationID: \"" + CompleteAmplificationID  + "\" }) return ss.UniqueID \n";
            System.out.println(command1);

            try (Session session = driver.session()) {
                StatementResult result = session.run(command1);
                ID = result.single().get(0).asString();

                System.out.println(CompleteAmplificationID+" -> "+ID);
            }

        } else{
            System.out.println("++++ ++++ I cannot determine the ID of "+file.getFileName().toString()+" ++++ ++++");
            System.exit(-1);
        }

        BufferedWriter logfile = null;
        try {
            logfile = new BufferedWriter(new FileWriter(ID+"_log.txt"));
            logfile.write("FILE "+file.getFileName().toString()+" "+ID+"\n\n");
        }  catch (IOException x) {
            System.err.println(x);
        }


        try {
            BufferedReader br = new BufferedReader(new FileReader(file.toString()));
            String line, lineaux;
            int num_insertion_file = 0;
            int numbrackets = 0;
            boolean isfirst = true;

            while ((line = br.readLine()) != null) {
                if (line.contains("ISCluster")) {
                    Insertion iscluster = new Insertion();
                    num_insertion_file++;
                    if (line.contains("{")) numbrackets++; //Dovrebbe essere sempre vero

                    while ((line = br.readLine()) != null) {
                        //Definisco il blocco
                        if (line.contains("{")) numbrackets++;
                        if (line.contains("}")) {
                            numbrackets--;
                            if (numbrackets == 0) {
                                //iscluster.PrintInsertion();
                                break;
                            }
                        }

                        if (line.contains("clu_id")) {
                            lineaux = line.substring(line.indexOf("\"") + 1);
                            lineaux = lineaux.substring(0, lineaux.indexOf("\""));
                            iscluster.clu_id = lineaux;
                        }

                        if (line.contains("cons_seq")) {
                            lineaux = line.substring(line.indexOf("\"") + 1);
                            if(lineaux.indexOf("\"") != -1) {
                                lineaux = lineaux.substring(0, lineaux.indexOf("\""));
                                iscluster.cons_seq = lineaux;
                            } else {
                                iscluster.cons_seq = lineaux;
                                do {
                                    lineaux = br.readLine();
                                    if(lineaux.indexOf("\"") != -1) iscluster.cons_seq = iscluster.cons_seq.concat(lineaux.substring(0, lineaux.indexOf("\"")));
                                    else iscluster.cons_seq = iscluster.cons_seq.concat(lineaux);
                                } while (lineaux.indexOf("\"") == -1);

                            }
                            //La fine è indicata da ". Potrebbe essere una linea sola o piu' linee
                            // di cui solo l'ultima contiene "
                            //System.out.println(num_insertion_file+" : "+iscluster.cons_seq);
                        }

                        if (line.contains("num_aln")) {
                            lineaux = line.substring(line.indexOf("l") + 3);
                            lineaux = lineaux.substring(0, lineaux.indexOf(","));
                            iscluster.num_aln = Integer.parseInt(lineaux);
                            //System.out.println(num_insertion_file+" : "+iscluster.num_aln);
                        }

                        if (line.contains("max_aln_score")) {
                            lineaux = line.substring(line.indexOf("e") + 2);
                            lineaux = lineaux.substring(0, lineaux.indexOf(","));
                            iscluster.max_aln_score = Integer.parseInt(lineaux);
                            //System.out.println(num_insertion_file+" : "+iscluster.max_aln_score);
                        }

                        if (line.contains("seq_len")) {
                            lineaux = line.substring(line.indexOf("n") + 2);
                            lineaux = lineaux.substring(0, lineaux.indexOf(","));
                            iscluster.seq_len = Integer.parseInt(lineaux);
                            //System.out.println(num_insertion_file+" : "+iscluster.seq_len);
                        }

                        if (line.contains("weight")) {
                            lineaux = line.substring(line.indexOf("t") + 2);
                            lineaux = lineaux.substring(0, lineaux.indexOf(","));
                            iscluster.weight = Integer.parseInt(lineaux);
                            //System.out.println(num_insertion_file+" : "+iscluster.weight);
                        }

                        if (line.contains("masked")) {
                            while (! (line = br.readLine()).contains("}") ) {
                                if (line.contains("\"")) {
                                    lineaux = line.substring(line.indexOf("\"") + 1);
                                    lineaux = lineaux.substring(0, lineaux.indexOf("\""));
                                    //System.out.println(num_insertion_file + " : " + lineaux);
                                    iscluster.masked.addElement(lineaux);
                                }
                                else {
                                    System.out.println("*** STRANA riga in masked : "+line);
                                    logfile.write("*** STRANA riga in masked : "+line+"\n");
                                    System.exit(-1);
                                }
                            }
                            numbrackets--;
                        }

                        if (line.contains("lab")) {
                            boolean stop = false;
                            String tid = null;
                            int c = 0, e = 0, o = 0, as = 0;
                            double cd = 0.0;
                            boolean isfloat = false;

                            while (!stop) {
                                line = br.readLine();

                                while (!line.contains("},") || line.contains("centroid")) {
                                    if (line.contains("target_id")) {
                                        lineaux = line.substring(line.indexOf("\"") + 1);
                                        lineaux = lineaux.substring(0, lineaux.indexOf("\""));
                                        tid = lineaux;
                                    }

                                    if (line.contains("centroid")) {
                                        lineaux = line.substring(line.indexOf("{") + 2);
                                        lineaux = lineaux.substring(0, lineaux.indexOf(","));
                                        try {
                                            c = Integer.parseInt(lineaux);
                                        } catch (NumberFormatException nfe) {
                                            cd = Double.parseDouble(lineaux);
                                            isfloat = true;
                                            System.out.println("*** ***\n\n"+line+" : "+cd);
                                        }
                                        //10^exp
                                        lineaux = line.substring(line.indexOf(",")+1);
                                        lineaux = lineaux.substring(lineaux.indexOf(",")+2);
                                        lineaux = lineaux.substring(0, lineaux.indexOf(" }"));
                                        e = Integer.parseInt(lineaux);
                                        if(isfloat) {
                                            cd *= Math.pow(10, e);
                                            c = Math.toIntExact(Math.round(cd));
                                            System.out.println(iscluster.clu_id+" : "+cd+" ^ "+e+" : "+c+" \n\n*** ***");
                                            isfloat = false;
                                        }
                                        else if (e > 0) {
                                            c *= Math.pow(10, e);
                                        }
                                        else if (e < 0) {
                                            cd = c * Math.pow(10, e);
                                            c = Math.toIntExact(Math.round(cd));
                                            System.out.println(iscluster.clu_id+" : "+cd+" ^ "+e+" : "+c+" \n\n*** ***");

                                        }
                                    }

                                    if (line.contains("offset")) {
                                        lineaux = line.substring(line.indexOf("t ") + 2);
                                        lineaux = lineaux.substring(0, lineaux.indexOf(","));
                                        o = Integer.parseInt(lineaux);
                                    }

                                    if (line.contains("aln_score")) {
                                        lineaux = line.substring(line.indexOf("e ") + 2);
                                        lineaux = lineaux.substring(0, lineaux.indexOf(","));
                                        as = Integer.parseInt(lineaux);
                                    }

                                    line = br.readLine();
                                    //System.out.println(line);
                                }
                                iscluster.targetlist.addElement(new Target(tid,c,o,as));
                                if(e < 0) System.out.println(iscluster.clu_id+" : "+tid+" "+c+" "+o+" "+as+" \n\n*** ***");

                                if ( (line = br.readLine()).contains("merged")) stop = true;
                                else if( !line.contains("{")) {
                                    System.out.println("*** STRANEZZE IN LAB "+ num_insertion_file + " : " + line);
                                    System.exit(-1);
                                }
                            }
                        }


                        //PROBLEMI
                        if (line.contains("ISCluster")) {
                            System.out.println("*** BLOCCO non chiuso");
                            logfile.write("*** BLOCCO non chiuso\n");
                            System.exit(-1);
                        }
                    }

                    //iscluster.PrintInsertion(ID);
                    if(isfirst) {
                        String command1 = "MERGE (ss:sample {UniqueID: \"" + ID + "\" })\n";

                        System.out.println("MATCH (ss:sample {UniqueID: \"" + ID + "\" }) return ss");
                        logfile.write("MATCH (ss:sample {UniqueID: \"" + ID + "\" }) return ss");
/*                        try (Session session = driver.session()) {
                            session.run(command1);
                        }
*/                        isfirst = false;
                    }

                    /***************************************************/
                    /***************************************************/
                    /***************************************************/
                    /***************************************************/
                    /***************************************************/

                    iscluster.AddInsertion(ID,driver,logfile,hsample);

                    /***************************************************/
                    /***************************************************/
                    /***************************************************/
                    /***************************************************/
                    /***************************************************/

                    try {
                        logfile.flush();
                    } catch (IOException e) {
                        e.printStackTrace();
                    }

                } // Iscluster
            }

            if(!hsample.isEmpty()) {

                Set<String> keys = hsample.keySet();
                Iterator<String> itr = keys.iterator();

                while (itr.hasNext()) {
                    String samplename = itr.next();
                    String command2 =  "MATCH (ss:sample {UniqueID: \"" +ID+ "\"})\n"+
                            "MATCH (ll:sample {UniqueID: \"" +samplename+ "\"})\n"+
                            "MERGE (ss)-[:gtris_sample {num: \""+hsample.get(samplename) +"\"}]-(ll)";

                    try (Session session = driver.session()) {
                        session.run(command2);
                        System.out.println("*** *** *** GTRIS between "+ID+ " and "+ samplename +"  ("+hsample.get(samplename)+" links)");
                        try {
                            logfile.write("*** *** *** GTRIS between "+ID+ " and "+ samplename+"  ("+hsample.get(samplename)+" links)\n");
                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                    }

                }

            }

            br.close();
        }
        catch (IOException x) {
            System.err.println(x);
        }

        driver.close();


        try {
            logfile.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

        System.out.println("DONE");
        System.exit(0);

    }
}
