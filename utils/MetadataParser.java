//https://neo4j.com/docs/driver-manual/1.7/

import org.neo4j.driver.v1.*;
import org.neo4j.driver.v1.exceptions.NoSuchRecordException;
import org.neo4j.driver.v1.util.Pair;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.lang.reflect.Field;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.StringTokenizer;

//import javafx.util.Pair;

class SubgSeq implements Comparable<SubgSeq>{
    String clu_id;
    int seq_len;
    String id;
    String cons_seq;
    int num_aln;


    public SubgSeq(String ci, int l, String d, String cs, int na) {
        clu_id = ci;
        seq_len = l;
        id = d;
        cons_seq = cs;
        num_aln = na;
    }

    public String toShortString() {
        return clu_id+" (len: "+seq_len+", aln: "+num_aln+")";
    }

    @Override
    public String toString() {
        return "[ " + seq_len + ", " + clu_id + ", " + cons_seq + ", " + num_aln +"]";
    }

    @Override
    public int compareTo(SubgSeq cmp) {
        return seq_len - cmp.seq_len;
    }
}


class CsvParser {
    String col_pos;
    String col_label;
    String new_label;
    String patient;
    String marker;
    String tissue;
    String timepoint;
    String ID;
}


public class MetadataParser {
    String ProjectID;
    String FUSIONID;
    String PoolID;
    String TagSequence;
    String SubjectID;
    String VectorType;
    String VectorID;
    String 	ExperimentID;
    String 	Tissue;
    String 	TimePoint;
    String 	DNAFragmentation;
    String 	PCRMethod;
    String 	TagIDextended;
    String 	Keywords;
    String 	CellMarker;
    String 	TagID;
    String 	NGSProvider;
    String 	NGSTechnology;
    String 	ConverrtedFilesDir;
    String 	ConverrtedFilesName;
    String 	SourceFileFolder;
    String 	SourceFileNameR1;
    String 	SourceFileNameR2;
    String 	DNAnumber;
    String 	ReplicateNumber;
    String 	DNAextractionDate;
    String 	DNAngUsed;
    String 	LinearPCRID;
    String 	LinearPCRDate;
    String 	SonicationDate;
    String 	LigationDate;
    String 	fstExpoPCRID; //era 1st
    String 	fstExpoPCRDate; //era 1st
    String 	sndExpoID; //era 2nd
    String 	sndExpoDate; //era 2st
    String 	FusionPrimerPCRID;
    String 	FusionPrimerPCRDate;
    String 	PoolDate;
    String 	SequencingDate;
    String 	VCN;
    String Genome;
    String 	SequencingRound;
    String 	Genotype;
    String 	TestGroup;
    String 	MOI;
    String Engraftment;
    String 	Transduction;
    String 	Notes;
    String 	AddedField1;
    String 	AddedField2;
    String 	AddedField3;
    String 	AddedField4;

    String concatenatePoolIDSeqRun; //Nuovo in 2020
    String  PoolID_SeqRun; //era PoolID-SeqRun
    String 	AddedField6_RelativeBloodPercentage;
    String 	AddedField7_PurityTestFeasibility;
    String 	AddedField8_FacsSeparationPurity;
    String 	Kapa;
    String 	ulForPool;
    String 	CompleteAmplificationID;
    String 	UniqueID;
    //Nuovi in 2020
    String 	StudyTestID;
    String 	StudyTestGroup;
    String 	MouseID;
    String 	Tigroup;
    String 	Tisource;
    String filename;
    //Nuovi per gli esperimenti sintetici
    String PathToFolderProjectID;
    String SamplesNameCheck;
    String TimepointDays;
    String TimepointMonths;
    String TimepointYears;
    //Nuovi 2022
    String ng_DNA_corrected;
    String Concatenate;
    String ISsoftware_name;
    String ISsoftware_version;
    String AdditionalISsoftware;



    void checkstring(Driver d) {

        String command1 = "match(n:subg) return n.clu_id, n.seq_len, n.cons_seq, n.num_aln";
        StatementResult result;
        ArrayList<SubgSeq> sequences = new ArrayList<SubgSeq>();


        String clu_id;
        int seq_len, num_aln;
        String id;
        String cons_seq;


        try (Session session = d.session()) {
            result = session.run(command1);

            while ( result.hasNext() ) {
                Record record = result.next();

                clu_id = record.get(0).asString();
                seq_len = Integer.parseInt(record.get(1).asString());
                cons_seq = record.get(2).asString();
                num_aln = Integer.parseInt(record.get(3).asString());
                id = clu_id.split("_")[0];

                sequences.add(new SubgSeq(clu_id,seq_len,id,cons_seq,num_aln));
            }


            Collections.sort(sequences);


            int sz = sequences.size();
            String candidate;
            int count;

            for(int i=0;i<sz;i++) {
                candidate = sequences.get(i).cons_seq;
                count = 0;

                for (int j = i + 1; j < sz; j++) {
                    if (sequences.get(j).cons_seq.contains(candidate)) {
                        count++;

                        command1 = "MATCH (a:subg {clu_id: \"" +sequences.get(j).clu_id+ "\"})-[l:gtris]-"+
                                "(b:subg {clu_id: \"" +sequences.get(i).clu_id+ "\"}) return l.case";

                        try (Session session2 = d.session()) {
                            result = session.run(command1);
                            if(result.hasNext())
                                System.out.println(sequences.get(i).toShortString()+ " -> "+ sequences.get(j).toShortString()+" GTRIS "+result.next().get(0).asString());
                            else System.out.println(sequences.get(i).toShortString()+ " -> "+sequences.get(j).toShortString());
                        }
                        //System.out.println(sequences.get(i).toString()+ " -> "+sequences.get(j).toString()+"\n");
                    }

                }
                if(count != 0) System.out.println("*** "+sequences.get(i).clu_id+ ": "+ count+"\n\n");
            }
        }



    }


    void insertN4j(Driver d) {

        int numvalues = 1;
        ///String command = "MERGE ("+this.UniqueID+":sample {";
        String command = "MERGE (n:sample {";

        command = command.replaceAll("-","_");

        for (Field field : this.getClass().getDeclaredFields()) {
            field.setAccessible(true);
            String name = field.getName();
            //System.out.println(name);
            Object value = null;
            try {
                value = field.get(this);
            } catch (IllegalAccessException e) {
                e.printStackTrace();
            }
            if(value != null){
                ///System.out.printf("%s: %s%n", name, value);
                if (!value.equals("empty") && !value.equals("NONE")) {
                    command = command.concat(name + ": \"" + value + "\", ");
                    numvalues++; //Info statistiche di valori significativi
                }
            }
        }

        command = command.substring(0,command.length()-2);
        command = command.concat("})");


        try (Session session = d.session())
        {
            // Wrapping Cypher in an explicit transaction provides atomicity
            // and makes handling errors much easier.
            try (Transaction tx = session.beginTransaction())
            {
                //tx.run("MERGE (a:exp {name: {x}})", parameters("x", name));
                tx.run(command);
                tx.success();  // Mark this write as successful.
                System.out.println("Inserted "+this.UniqueID);
            }
        }
    }


    static MetadataParser getN4j(Driver d, CsvParser csvp) {
        try (Session session = d.session())
        {

            MetadataParser mp = null;

            Field[] f = MetadataParser.class.getDeclaredFields();

            String returnres = "return ";

            for(int i=0;i<f.length;i++) returnres+= "n."+f[i].getName()+", ";

            returnres = returnres.substring(0,returnres.length()-2);

            StatementResult result = session.run("MATCH (n:sample {UniqueID: \""+csvp.ID+"\" }) "+returnres);

            if (result.hasNext()) { //forse anche != null
                mp = new MetadataParser();
                Record res = null;

                try {
                    res = result.single();
                } catch (NoSuchRecordException e) {
                    System.out.println("******* "+csvp.ID+" *******");
                    e.printStackTrace();
                }

                List<Pair<String,Value>> values = res.fields();
                int j = 0;
                for (Pair<String, Value> pair : values) {

                    try {
                        f[j++].set(mp,pair.value().asString());
                    } catch (IllegalAccessException e) {
                        e.printStackTrace();
                    }
                    ///System.out.println(pair.key()+" -- "+pair.value());

                }

                //A volte UTR VS UT */
                if(!csvp.patient.equalsIgnoreCase(mp.SubjectID) || /*!csvp.tissue.equalsIgnoreCase(mp.Tissue) ||*/ (Integer.parseInt(csvp.timepoint) != Integer.parseInt(mp.TimePoint)) ) {
                    System.out.println("*** DATI DIFFORMI " + csvp.ID + " " + csvp.patient + " " + mp.SubjectID + " " + csvp.tissue + " " + mp.Tissue+ " " + Integer.parseInt(csvp.timepoint) + " " + Integer.parseInt(mp.TimePoint));
                    mp = null;
                } 
                else if( !mp.filename.equalsIgnoreCase("tobereplaced")) {//Ho gia' un filename
                    if(!mp.filename.equalsIgnoreCase(csvp.col_label) ) {
                        System.out.println("*** FILENAME ESISTENTE E DIVERSO " + csvp.ID + " :-"+mp.filename+ "-:-"+ csvp.col_label +"-");
                        mp = null;
                    }
                }

            } else {
                System.out.println("*** Non esiste " + csvp.ID + " !");
            }


            return(mp);
        }

    }

    static void extendN4j(Driver d, MetadataParser mp, String filename) {
        try (Session session = d.session())
        {
            String query = "MATCH (n:sample {UniqueID: \""+mp.UniqueID+"\" }) ";
            query += "SET n.filename = \""+filename+"\" return n";
            StatementResult result = session.run(query);

            if (result != null) {
                System.out.println("Aggiunto per "+mp.UniqueID);
            } else System.out.println("ERRORE con query; "+query);

        }
    }


    public static void main(String[] args) {

        //Command LINE
        Path file = null;

        String[] filename;
        MetadataParser mp;
        Driver driver;


        driver = GraphDatabase.driver("bolt://localhost:7687", AuthTokens.basic(/*UID*/, /*PWD*/));

        if (file == null)
            if (args.length == 0) {
                System.out.println("Indicare il file da leggere! \n Uso: java MetadataParser <file with a complete path> | checkstring");
                System.exit(-1);
            } else if(args[0].compareTo("checkstring") == 0) {

                mp = new MetadataParser();
                mp.checkstring(driver);
                driver.close();
                System.exit(0);

            } else file = Paths.get(args[0]);



        System.out.println("Sto leggendo "+file);
        String fext = file.toString();
        fext = fext.substring(fext.length()-3);


        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if ( fext.compareToIgnoreCase ("tsv") == 0) {

            try {


                BufferedReader br = new BufferedReader(new FileReader(file.toString()));
                String line;
                StringTokenizer colname = null, lineaux;
                Field[] f = null;
                int j,k = 0,numelem=0;

                while ((line = br.readLine()) != null) {
                    //System.out.println(line);
                    mp = new MetadataParser();

                    /*** *** *** *** Riga intestazione delle colonne *** *** *** ***/
                    if (line.contains("1stExpoPCRDate")) {
                        line = line.replaceAll("1st","fst");
                        line = line.replaceAll("2nd","snd");
                        line = line.replaceAll("concatenate PoolID-SeqRun","PoolID_SeqRun");
                        colname = new StringTokenizer(line," \t"); //Intestazione delle colonne.

                        k = colname.countTokens();
                        f = new Field[k];
                        System.out.println("Ho "+k+" colonne");
                        for(int i=0;i<k;i++) {

                            try {
                                f[i] = MetadataParser.class.getDeclaredField(colname.nextToken());
                                ///System.out.println(f[i].getName());
                            } catch (NoSuchFieldException e) {
                                e.printStackTrace();
                            }
                        }
                        ///System.out.println("\n *** *** *** *** *** \n");
                        continue;
                    } //Fine riga intestazione

                    numelem++;
                    line = line.replaceAll("\t\t","\tempty\t"); //Ho serie di \t che senno' non prende correttamente
                    line = line.replaceAll("\t\t","\tempty\t");

                    lineaux = new StringTokenizer(line,"\t");
                    ///System.out.println("Ho "+lineaux.countTokens()+" colonne");
                    j = 0;
                    if (lineaux.countTokens() < k-1) System.err.println(line+"\n*** Ho solo "+lineaux.countTokens()+" colonne invece di "+(k-1));
                    if (lineaux.countTokens() > k-1) System.err.println(line+"\n*** Ho ben "+lineaux.countTokens()+" colonne invece di "+(k-1));
                    while ( lineaux.hasMoreTokens() ) {
                        try {
                            f[j++].set(mp,lineaux.nextToken());
                            ///String a = lineaux.nextToken();
                            ///System.out.println("Inserisco "+a+" in "+f[j]);
                            ///f[j++].set(mp,a);
                            ///System.out.println(f[j-1].getName()+" : "+f[j-1].get(mp));
                        } catch (IllegalAccessException e) {
                            e.printStackTrace();
                        }

                    }
                    mp.filename = "tobereplaced";

                    ///System.out.print("Row "+numelem+ " ");
                    mp.insertN4j(driver);

                    ///System.out.println("\n *** *** *** *** *** \n");
                    //mp.printMP();
                    //System.exit(0);
                }
                br.close();

                //mp.printN4j(driver);
                driver.close();
            } catch (IOException x) {
                System.err.println(x);
            } } //FINE FILE TSV

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        else if (fext.compareToIgnoreCase("csv") == 0) {
            try {
                BufferedReader br = new BufferedReader(new FileReader(file.toString()));
                String line;
                StringTokenizer colname = null, lineaux;
                Field[] f = null;
                int j,k = 0;

                CsvParser csvp;

                //Leggo le colonne
                //NOTA: nella conversione Excel - csv ha inserito un carattere strano all'inizio. Va eliminato con nano
                line = br.readLine();
                //System.out.println(line);

                lineaux = new StringTokenizer(line,";");

                ///Se la riga inizia con ; dovrei partire da 1 con il for i
                k = lineaux.countTokens();
                f = new Field[k];
                ///System.out.println("Ho "+k+" colonne");

                /*** *** *** *** Riga intestazione delle colonne *** *** *** ***/
                for(int i=0;i<k;i++) {
                    try {
                        f[i] = CsvParser.class.getDeclaredField(lineaux.nextToken());
                        //System.out.println(f[i].getName());
                    } catch (NoSuchFieldException e) {
                        e.printStackTrace();
                    }
                }

                while ((line = br.readLine()) != null) {
                    ///System.out.println(line);
                    csvp = new CsvParser();
                    lineaux = new StringTokenizer(line,",;");

                    j = 0;
                    mp = null;
                    while ( lineaux.hasMoreTokens() ) {
                        try {
                            f[j++].set(csvp,lineaux.nextToken());
                            ///String a = lineaux.nextToken();
                            ///System.out.println("Inserisco "+a+" in "+f[j]);
                            ///f[j++].set(mp,a);
                            ///System.out.println(f[j-1].getName()+" : "+f[j-1].get(csvp));
                        } catch (IllegalAccessException e) {
                            e.printStackTrace();
                        }
                    }

                    mp = MetadataParser.getN4j(driver,csvp);

                    if (mp == null)  continue;
                    MetadataParser.extendN4j(driver,mp,csvp.col_label);

/////////////////////////////////////////////////
/*
                    lineaux = new StringTokenizer(line,",;");

                    j = 0;
                    if(k > lineaux.countTokens()) {
                        j++;
                        mp.col_pos = "-1";
                    }

                    //while (lineaux.hasMoreTokens()) System.out.print(lineaux.nextToken()+" -- ");
                    //System.out.println();
                    while (lineaux.hasMoreTokens()) {
                            try {
                                f[j++].set(mp,lineaux.nextToken());
                                //System.out.println(f[j-1].getName()+" : "+f[j-1].get(mp));
                            } catch (IllegalAccessException e) {
                                e.printStackTrace();
                            }

                    }

                    //mp.getN4j(driver,mp.ID);
                    mp.printN4j(driver,mp.ID);
                    //mp.extendN4j(driver);
                    mp.printMP();

                    System.exit(0);
*/
                }

            } catch (IOException x) {
                System.err.println(x);
            } }
        else {
            System.out.println("Indicare un file .csv o .tsv, ho "+file.toString()+" "+fext);
            System.exit(-2);
        }



        System.out.println("Fatto");
        System.exit(0);
    }

}


/*

| (:sample {VCN: "1.98", UniqueID: "ID00000000000000003304", SequencingDate: "2012-01-18", Keywords: "CD34", PoolID: "MLD01-POOL12-GSK", LinearPCRDate: "2012-01-24", ProjectID: "MLD", AddedField7_PurityTestFeasibility: "1", DNAngUsed: "100", DNAextractionDate: "2008-01-01", ConverrtedFilesName: "TTGCCG.fa.gz", FUSIONID: "FB#334.32", SourceFileNameR1: "MLD01-12_S2_L001_R1_001.fastq.gz", Tissue: "BM", SourceFileNameR2: "MLD01-12_S2_L001_R2_001.fastq.gz", fstExpoPCRID: "MM#35.10", CellMarker: "CD34", VectorType: "lenti", DNAnumber: "MLD01-DNA205", TagID: "TTGCCG", ConverrtedFilesDir: "GSK-IlluminaTo454/MLD01-pool12-GSK", CompleteAmplificationID: "MLD_MLD01-POOL12-GSK_TTGCCG_MLD01_MLD01-DNA205_lenti_751.pCCLsin.PPT.hPGK.ARSA.Wpre-mut-KANA_NONE_BM_1_NONE_LAM-PCR_0540_CD34_CD34_GSK-IlluminaTo454_NONE_NONE_NONE", fstExpoPCRDate: "2012-01-25", DNAFragmentation: "ACI", PoolDate: "2012-05-29", TimePoint: "0540", TagSequence: "TTGCCG", TagIDextended: "TTGCCG", sndExpoID: "MM#36.10", AddedField6_RelativeBloodPercentage: "1.77", SourceFileFolder: "ClinicalTrials/Projects/MLD/data/raw/20121018_GSK_Illumina/MLD01_POOL12/MLD01_POOL12", FusionPrimerPCRID: "FB#334.32", FusionPrimerPCRDate: "2012-03-23", NGSProvider: "GSK-IlluminaTo454", PoolID_SeqRun: "MLD01-POOL12-GSK-1", NGSTechnology: "MiSeq", sndExpoDate: "2012-01-25", ReplicateNumber: "1", LinearPCRID: "MM#34.10", PCRMethod: "LAM-PCR", SubjectID: "MLD01", Genome: "HUMAN", SequencingRound: "1", VectorID: "751.pCCLsin.PPT.hPGK.ARSA.Wpre-mut-KANA"}) 
 */