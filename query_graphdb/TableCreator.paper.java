import org.neo4j.driver.v1.*;
import org.neo4j.driver.v1.exceptions.value.LossyCoercion;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.text.SimpleDateFormat;
import java.lang.Math.*;


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


public class TableCreator {


    //Linko le subregion che sono sottostringhe di altre
    static void checkstring(Driver d) {

        BufferedWriter logfile = null;
        try {
            logfile = new BufferedWriter(new FileWriter(new SimpleDateFormat("'Substring-'yyyyMMddHHmm'.txt'").format(new Date()) ));
        }  catch (IOException x) {
            System.err.println(x);
        }



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
                //seq_len = Integer.parseInt(record.get(1).asString());
                seq_len = record.get(1).asInt();
                cons_seq = record.get(2).asString();
                //num_aln = Integer.parseInt(record.get(3).asString());
                num_aln = record.get(3).asInt();
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
                            if(result.hasNext()) {
                                String s = sequences.get(i).toShortString() + " -> " + sequences.get(j).toShortString() + " GTRIS " + result.next().get(0).asString();
                                System.out.println(s);
                                logfile.write(s+"\n");
                            }
                            else {
                                String s = sequences.get(i).toShortString()+ " -> "+sequences.get(j).toShortString();
                                System.out.println(s);
                                logfile.write(s+"\n");
                            }
                        } catch (IOException e) {
                            e.printStackTrace();
                        }

                    //System.out.println(sequences.get(i).toString()+ " -> "+sequences.get(j).toString()+"\n");
                    }

                }

                try {
                    if(count != 0) {
                        String s = "*** " + sequences.get(i).clu_id + ": " + count + "\n\n";
                        System.out.println(s);
                        logfile.write(s);
                        logfile.flush();
                    }
                } catch (IOException e) {
                    e.printStackTrace();
                }

            }
        }


        try {
            logfile.flush();
            logfile.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


/************************** ************************** ************************** ************************** **************************
 ************************** ************************** ************************** ************************** **************************
 ************************** ************************** ************************** ************************** **************************
 ************************** ************************** ************************** ************************** **************************
 ************************** ************************** ************************** ************************** **************************
 ************************** ************************** ************************** ************************** **************************/


    static void table(Driver d, float threshold, String wlabel) {

        BufferedWriter logfile = null;
        try {
            //logfile = new BufferedWriter(new FileWriter(new SimpleDateFormat("'Link-'yyyyMMddHHmm'.txt'").format(new Date()) ));
            logfile = new BufferedWriter(new FileWriter("TableFile.log"));
        }  catch (IOException x) {
            System.err.println(x);
        }

        BufferedWriter tablefile_u = null;
        BufferedWriter tablefile_s = null;
        BufferedWriter tablefile_t = null;
        BufferedWriter tablefile_c = null;


        if(wlabel.compareToIgnoreCase("all") == 0) {
            try {
                tablefile_u = new BufferedWriter(new FileWriter("TableFile_clu.csv" ));
                tablefile_s = new BufferedWriter(new FileWriter("TableFile_shear.csv" ));
                tablefile_t = new BufferedWriter(new FileWriter("TableFile_tag.csv" ));
                tablefile_c = new BufferedWriter(new FileWriter("TableFile_combo.csv" ));
            }  catch (IOException x) {
                System.err.println(x);
            }


        } else {
            String ext = wlabel;
            if(wlabel.compareToIgnoreCase("weight") == 0) ext = "clu";

            try {
                tablefile_u = new BufferedWriter(new FileWriter("TableFile_"+ext+".csv" ));
            }  catch (IOException x) {
                System.err.println(x);
            }
        }


        ArrayList<String> subg = new ArrayList<String>();


        //https://www.geeksforgeeks.org/arraylist-contains-java/
        Hashtable<String, Integer> htwgt_u = new Hashtable<String, Integer>();
        Hashtable<String, Integer> htwgt_s = new Hashtable<String, Integer>();
        Hashtable<String, Integer> htwgt_t = new Hashtable<String, Integer>();
        Hashtable<String, Integer> htwgt_c = new Hashtable<String, Integer>();

        Hashtable<String, String> genpos = new Hashtable<String, String>();

        //Per prima cosa determino il numero massimo di righe -> subg / integrazioni NO, non serve
        //Poi il numero massimo di colonne -> ID / esperimenti / campioni

        String commandcols = "match (n:sample) return count(n)";
        int numcols = 0;

        String commandnamecols = "match (n:sample) return n.UniqueID order by n.UniqueID";


        ArrayList<String> namecols = new ArrayList<String>();

        StatementResult result;
        //I Create the resulting matrix
        try (Session session = d.session()) {

            result = session.run(commandcols);
            numcols = result.single().get(0).asInt();

            result = session.run(commandnamecols);
            while (result.hasNext()) {
                Record record = result.next();
                namecols.add(record.get(0).asString());
            }
        }

        System.out.println("I have "+numcols+" experiments, threshold "+threshold);
        try {
            logfile.write("I have "+numcols+" experiments, threshold "+threshold+"\n");
        } catch (IOException e) {
            e.printStackTrace();
        }
        int[] tablefilerow_u = new int[numcols];
        int[] tablefilerow_s = new int[numcols];
        int[] tablefilerow_t = new int[numcols];
        int[] tablefilerow_c = new int[numcols];



        //Header of the table
        //Mancano,Chr,Pos,Strand,
        try {
            tablefile_u.write("ID,Chr,Pos,Strand,numsubg,numexp,SUBG,");
            for (int j = 0; j < numcols-1; j++) {
                tablefile_u.write(namecols.get(j) + ",");
            }
            tablefile_u.write(namecols.get(numcols-1) +"\n");
            tablefile_u.flush();

            //I write 4 files
            if(wlabel.compareToIgnoreCase("all") == 0) {
                tablefile_s.write("ID,Chr,Pos,Strand,numsubg,numexp,SUBG,");
                for (int j = 0; j < numcols-1; j++) {
                    tablefile_s.write(namecols.get(j) + ",");
                }
                tablefile_s.write(namecols.get(numcols-1) +"\n");
                tablefile_s.flush();

                tablefile_t.write("ID,Chr,Pos,Strand,numsubg,numexp,SUBG,");
                for (int j = 0; j < numcols-1; j++) {
                    tablefile_t.write(namecols.get(j) + ",");
                }
                tablefile_t.write(namecols.get(numcols-1) +"\n");
                tablefile_t.flush();

                tablefile_c.write("ID,Chr,Pos,Strand,numsubg,numexp,SUBG,");
                for (int j = 0; j < numcols-1; j++) {
                    tablefile_c.write(namecols.get(j) + ",");
                }
                tablefile_c.write(namecols.get(numcols-1) +"\n");
                tablefile_c.flush();
            }



        } catch (IOException e) {
            e.printStackTrace();
        }


        String commandsubg;

        if(wlabel.compareToIgnoreCase("all") == 0) commandsubg = "match(n:subg) return n.clu_id, n.weight, n.shear, n.tag, n.combo order by n.clu_id";
        else commandsubg = "match(n:subg) return n.clu_id, n."+wlabel+" order by n.clu_id";

        int numrows = 0;

        ArrayList<String> namerows = new ArrayList<String>();


        String qres;

        /*
         *
         *
         * CICLO per collezionare le subg
         *
         *
         */
        try (Session session = d.session()) {
            result = session.run(commandsubg);

            while (result.hasNext()) {
                Record record = result.next();
                String res0 = record.get(0).asString();
                subg.add(res0);
                //System.out.println(record.get(0).asString()+" "+record.get(1).asString());

                int val;

                if(wlabel.compareToIgnoreCase("all") == 0) {

                    htwgt_u.put(res0, record.get(1).asInt());
                    htwgt_s.put(res0, record.get(2).asInt());
                    htwgt_t.put(res0, record.get(3).asInt());
                    htwgt_c.put(res0, record.get(4).asInt());

                } else {
                    htwgt_u.put(res0, record.get(1).asInt());
                }
            }
        }

        Collections.sort(subg);

        ArrayList<String> nameofrows = new ArrayList<String>();
        ArrayList<String> htris = new ArrayList<String>();
        ArrayList<String> htris1 = new ArrayList<String>();
        ArrayList<String> htris2 = new ArrayList<String>();



        String candidate1 ="",candidate2;
        String prefixclu1, prefixclu2;
        String commandgtris, commandtarget;
        int count, diameter;
        int numsubg;
        //, numprobl;

        //Inizializzo a casa senno' scoccia
        String target ="T", offset="O";
        int seqlen = 0, centroid=-1, numaln;

        /* Potenziali candidati */
        String ptarget, poffset;
        int pcentroid;
        boolean newtarget = false;


        /*
         *
         *
         * CICLO per raggruppare le subg
         *
         *
         */

        //for(int i=0;i<subg.size();i++) {
        while(subg.size() > 0) {

            System.out.println(subg.size() + " " + subg.get(0));

            try {
                logfile.write("\n\n"+subg.size() + " " + subg.get(0)+"\n");
            } catch (IOException e) {
                e.printStackTrace();
            }

            candidate1 = subg.get(0);
            prefixclu1 = candidate1.split("_")[0];

            nameofrows.add(candidate1);
            numrows++;

            Arrays.fill(tablefilerow_u, 0);
            Arrays.fill(tablefilerow_s, 0);
            Arrays.fill(tablefilerow_t, 0);
            Arrays.fill(tablefilerow_c, 0);

            htris.clear();
            htris1.clear();
            htris2.clear();


            try (Session session = d.session()) {

                /************************** Inizializzo il target **************************/
                commandtarget = "match(n:subg {clu_id: '"+candidate1+"' })  return n.num_aln, n.seq_len";
                //System.out.println(commandtarget);

                result = session.run(commandtarget);
                Record r = result.next();

                numaln = r.get(0).asInt();
                seqlen = r.get(1).asInt();


                if( numaln == 1)  {
                    commandtarget = "match(n:subg {clu_id: '"+candidate1+"' })-[l]-(m:pos)  return m.target_id, m.centroid, l.offset ";
                    //System.out.println(commandtarget);

                    result = session.run(commandtarget);
                    r = result.next();

                    target = r.get(0).asString();
                    centroid = r.get(1).asInt();
                    if (r.get(2).asInt() >= 0) offset = "+";
                    else offset = "-";

                    logfile.write(target+" "+centroid+" "+offset+" ("+numaln+","+seqlen+")\n");

                } else  {
                    //Nuove query
                    //match (s:subg)-[l:main_insertion]-(p:pos) where l.val < 0.4 return p.target_id, p.centroid, l.offset, s.clu_id;

                    //match(n:subg {clu_id: 'ID00000000000000006244_subg_4' })-[l:gtris]-(m) where l.case <> "4" return m.clu_id, l.case, m.num_aln, m.seq_len, l.val order by l.case, m.clu_id;
                    //match(n:subg {clu_id: 'ID00000000000000006244_subg_4' })-[l:gtris]-(m) where l.case = "4" and l.val < 0.3 return m.clu_id, l.case, m.num_aln, m.seq_len, l.val order by l.case, m.clu_id;

                    commandtarget = "match(s:subg {clu_id: '"+candidate1+"'})-[l:main_insertion]-(p:pos) where l.val < "+threshold+" return  p.target_id, p.centroid, l.offset, l.val";
                    //ERA
                    //commandtarget = "match(n:subg {clu_id: '"+candidate1+"'})-[l]-(m:pos)  return l.aln_score, m.target_id, m.centroid, l.offset  order by l.aln_score desc limit 2";
                    result = session.run(commandtarget);

                    float aa;

                    if(result.hasNext()) {
                        r = result.next();

                        target = r.get(0).asString();
                        centroid = r.get(1).asInt();
                        if (r.get(2).asInt() >= 0) offset = "+"; else offset = "-";
                        //aa= r.get(3).asFloat();
                        try {
                            aa= r.get(3).asFloat();
                        } catch (LossyCoercion e) {
                            aa= (float) r.get(3).asDouble();
                        }

                    } else {
                        target = "R";
                        centroid = -1;
                        offset = "R";
                        aa = -1; //Useless, just for logfile
                    }

                    logfile.write(target+" "+centroid+" "+offset+" ("+numaln+","+seqlen+"): "+aa+" \n");

                }




                /************************** Ciclo sul suo primo intorno **************************/
                commandgtris = "match(n:subg {clu_id: '" + candidate1 + "' })-[l:gtris]-(m) return m.clu_id, l.case, m.num_aln, m.seq_len, l.val order by m.clu_id";
                result = session.run(commandgtris);

                while (result.hasNext()) {
                    Record record = result.next();

                    String gcase = record.get(1).asString();
                    if(gcase.compareTo("4") == 0 )
                        try {
                            if (record.get(4).asFloat() >= threshold)
                                continue; //l.val is not NULL only for GTRIS 4
                        } catch (LossyCoercion e) {
                            if (record.get(4).asDouble() >= threshold)
                                continue; //l.val is not NULL only for GTRIS 4
                        }

                    //From here on GTRIS 1-3 and 4 but < threshold
                    htris.add(record.get(0).asString());
                    logfile.write(record.get(0).asString()+"  (gt "+record.get(1).asString()+") aln: "+record.get(2).asInt()+" seqlen: "+record.get(3).asInt()+"\n");

                    if(record.get(3).asInt() > seqlen  || ( record.get(3).asInt() == seqlen && record.get(2).asInt() > 1)  ) {
                        seqlen = record.get(3).asInt();

                        //num_aln - non so se succedera' mai - Invece si
                        if(record.get(2).asInt() == 1) {
                            commandtarget = "match(n:subg {clu_id: '" + record.get(0).asString() + "' })-[l]-(m:pos)  return m.target_id, m.centroid, l.offset ";
                            StatementResult result1 = session.run(commandtarget);
                            Record r1 = result1.next();

                            ptarget = r1.get(0).asString();
                            pcentroid = r1.get(1).asInt();
                            if (r1.get(2).asInt() >= 0) poffset = "+";
                            else poffset = "-";

                            if(target.equalsIgnoreCase(ptarget) && offset.equalsIgnoreCase(poffset) && java.lang.Math.abs(centroid-pcentroid) < 3) newtarget = false;
                            else {
                                newtarget = true;
                                target = ptarget;
                                centroid = pcentroid;
                                offset = poffset;
                                logfile.write("Very strange\n");
                            }

                        } else {

                            //If there is no main target it's a repeat.
                            commandtarget = "match(s:subg {clu_id: '"+record.get(0).asString()+"'})-[l:main_insertion]-(p:pos) where l.val < "+threshold+" return  p.target_id, p.centroid, l.offset";

                            //commandtarget = "match(n:subg {clu_id: '"+record.get(0).asString()+"'})-[l]-(m:pos)  return m.target_id, m.centroid, l.offset  order by l.aln_score desc limit 1";
                            StatementResult result1 = session.run(commandtarget);

                            if(result1.hasNext()) {
                                Record r1 = result1.next();
                                ptarget = r1.get(0).asString();
                                pcentroid = r1.get(1).asInt();
                                if (r1.get(2).asInt() >= 0) poffset = "+"; else poffset = "-";

                                if(target.equalsIgnoreCase(ptarget) && offset.equalsIgnoreCase(poffset) && java.lang.Math.abs(centroid-pcentroid) < 3) newtarget = false;
                                else {
                                    if(!target.equalsIgnoreCase("R")) newtarget = true;
                                    target = "R";
                                    centroid = -1;
                                    offset = "R";
                                }
                            } else {
                                if(!target.equalsIgnoreCase("R")) newtarget = true;
                                target = "R";
                                centroid = -1;
                                offset = "R";
                            }

                            /*
                            Record r1 = result1.next();

                            ptarget = r1.get(0).asString();
                            pcentroid = r1.get(1).asInt();
                            if (r1.get(2).asInt() >= 0) poffset = "+";
                            else poffset = "-";

                            if(target.equalsIgnoreCase(ptarget) && offset.equalsIgnoreCase(poffset) && java.lang.Math.abs(centroid-pcentroid) < 3) newtarget = false;
                            else {
                                newtarget = true;
                                target = "R";
                                centroid = -1;
                                offset = "R";
                            }
                            */
                        }
                        if(newtarget) {
                            logfile.write("TARGET CHANGED in "+target+" "+centroid+" "+offset+" \n");
                            newtarget = false;
                        }
                    }
                }
            } catch (IOException e) {
                e.printStackTrace();
            }


            //Inizializzo l'array per cercare ulteriori relazioni (ovvero diametro > 1)

            htris1.addAll(htris);
            htris.add(candidate1); //Sennò non funziona htris.contains(candidate2) per le cricche

            System.out.println("Start with " + htris.size()+ " subg");
            numsubg = htris.size();

            ///System.out.flush();
            try {
                logfile.write(candidate1+" start with: "+htris.size()+"\n");
            } catch (IOException e) {
                e.printStackTrace();
            }


            diameter = 0;


            /************************** Ciclo sugli eventuali successivi intorni **************************/
            do {
                count = 0;
                htris2.clear();

                for (int j = 0; j < htris1.size(); j++) {
                    //commandgtris = "match(n:subg {clu_id: '" + htris1.get(j) + "' })-[l:gtris]-(m) return m.clu_id, l.case, m.num_aln, m.seq_len order by m.clu_id";
                    commandgtris = "match(n:subg {clu_id: '" + htris1.get(j) + "' })-[l:gtris]-(m) return m.clu_id, l.case, m.num_aln, m.seq_len, l.val order by m.clu_id";


                    try (Session session = d.session()) {

                        result = session.run(commandgtris);

                        while (result.hasNext()) {
                            Record record = result.next();

                            String gcase = record.get(1).asString();
                            if(gcase.compareTo("4") == 0 )
                                try {
                                    if (record.get(4).asFloat() >= threshold)
                                        continue; //l.val is not NULL only for GTRIS 4
                                } catch (LossyCoercion e) {
                                    if (record.get(4).asDouble() >= threshold)
                                        continue; //l.val is not NULL only for GTRIS 4
                                }
                            candidate2 = record.get(0).asString();

                            //Nuova subg, quindi il diametro cresce. Controllo in tutti e 3 i set, a volte trovo più volte lo stesso...
                            if (!htris.contains(candidate2) && !htris1.contains(candidate2) && !htris2.contains(candidate2)) {
                                ///System.out.println("Aggiungo: "+candidate2);

                                count++;
                                htris2.add(candidate2);
                                logfile.write(record.get(0).asString()+"  (gt "+record.get(1).asString()+") aln: "+record.get(2).asInt()+" seqlen: "+record.get(3).asInt()+"\n");

                                /*if(Integer.parseInt(record.get(3).asString()) > seqlen  || ( Integer.parseInt(record.get(3).asString()) == seqlen && Integer.parseInt(record.get(2).asString()) > 1)  ) {
                                    seqlen = Integer.parseInt(record.get(3).asString());*/
                                if(record.get(3).asInt() > seqlen  || ( record.get(3).asInt() == seqlen && record.get(2).asInt() > 1)  ) {
                                    seqlen = record.get(3).asInt();

                                    if(record.get(2).asInt() == 1) {
                                        commandtarget = "match(n:subg {clu_id: '" + candidate2 + "' })-[l]-(m:pos)  return m.target_id, m.centroid, l.offset ";
                                        StatementResult result2 = session.run(commandtarget);
                                        Record r2 = result2.next();

                                        ptarget = r2.get(0).asString();
                                        pcentroid = r2.get(1).asInt();
                                        if (r2.get(2).asInt() >= 0) poffset = "+";
                                        else poffset = "-";

                                        if(target.equalsIgnoreCase(ptarget) && offset.equalsIgnoreCase(poffset) && java.lang.Math.abs(centroid-pcentroid) < 3) newtarget = false;
                                        else {
                                            newtarget = true;
                                            target = ptarget;
                                            centroid = pcentroid;
                                            offset = poffset;
                                            logfile.write("Very very strange\n");
                                        }

                                    } else {
                                        //commandtarget = "match(n:subg {clu_id: '"+record.get(0).asString()+"'})-[l]-(m:pos)  return m.target_id, m.centroid, l.offset  order by l.aln_score desc limit 1";
                                        commandtarget = "match(s:subg {clu_id: '"+record.get(0).asString()+"'})-[l:main_insertion]-(p:pos) where l.val < "+threshold+" return  p.target_id, p.centroid, l.offset";

                                        StatementResult result2 = session.run(commandtarget);

                                        if(result2.hasNext()) {
                                            Record r2 = result2.next();
                                            ptarget = r2.get(0).asString();
                                            pcentroid = r2.get(1).asInt();
                                            if (r2.get(2).asInt() >= 0) poffset = "+";
                                            else poffset = "-";

                                            if (target.equalsIgnoreCase(ptarget) && offset.equalsIgnoreCase(poffset) && java.lang.Math.abs(centroid - pcentroid) < 3)
                                                newtarget = false;
                                            else {
                                                if (!target.equalsIgnoreCase("R")) newtarget = true;
                                                target = "R";
                                                centroid = -1;
                                                offset = "R";
                                            }
                                        } else {
                                                if(!target.equalsIgnoreCase("R")) newtarget = true;
                                                target = "R";
                                                centroid = -1;
                                                offset = "R";
                                            }
                                    }

                                    if(newtarget) {
                                        logfile.write("TTARGET CHANGED in "+target+" "+centroid+" "+offset+" \n");
                                        newtarget = false;
                                    }
                                }

                            }

                        }
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }

                if (count != 0) {
                    //La prima volta non va fatto
                    if (diameter != 0) htris.addAll(htris1);
                    htris1.clear();
                    htris1.addAll(htris2); //Considero solo le nuove subg
                    diameter++;

                    System.out.println("New subg: " + htris2.size()+" (diam: "+diameter+")");
                    numsubg += htris2.size();
                    try {
                        logfile.write("New sets: " + htris2.size()+" "+diameter+"\n");
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
            } while (count != 0);

            if (diameter != 0) {
                htris.addAll(htris1);
                System.out.println("End with " + htris.size()+ " subg, diameter "+diameter);
            } else System.out.println("End with " + htris.size());

            //if (numprobl!=0) System.out.println("PROBLEM: "+numprobl);


            /*
             *
             *
             * NEW ROW
             *
             *
             */

            int not0 = 0;

            for (int j = 0; j < htris.size(); j++) {
                candidate2 = htris.get(j);
                prefixclu2 = candidate2.split("_")[0];

                //System.out.println(prefixclu2);

                try {
                    tablefilerow_u[namecols.indexOf(prefixclu2)] += htwgt_u.get(candidate2);
                    tablefilerow_s[namecols.indexOf(prefixclu2)] += htwgt_s.get(candidate2);
                    tablefilerow_t[namecols.indexOf(prefixclu2)] += htwgt_t.get(candidate2);
                    tablefilerow_c[namecols.indexOf(prefixclu2)] += htwgt_c.get(candidate2);
                } catch(NullPointerException e) {
                    System.out.println("NULL "+namecols.indexOf(prefixclu2)+" "+candidate2+" "+htwgt_u.contains(candidate2));
                    try {
                        logfile.flush();
                    } catch (IOException e1) {
                        e1.printStackTrace();
                    }
                    System.exit(-1);
                }

                subg.remove(candidate2);
            }


            try {
                logfile.write((numrows-1)+","+numsubg+","+nameofrows.get(numrows - 1) + "\n");

                logfile.write("FINAL: "+target+" "+centroid+" "+offset+" Diam: "+diameter);

                String poslabel = target+" "+centroid;
                if(genpos.contains(poslabel)) logfile.write("*** Duplicated entry: "+target+" "+centroid+" : "+genpos.get(poslabel)+"\n");
                else genpos.put(poslabel,candidate1);


               // "ID,Chr,Pos,Strand,numsubg,numexp,SUBG,"
                for (int j = 0; j < numcols; j++) if(tablefilerow_u[j] != 0) not0++;

                tablefile_u.write(numrows-1+","+target+","+centroid+","+offset+","+numsubg+","+not0+","+nameofrows.get(numrows - 1) + ",");
                tablefile_s.write(numrows-1+","+target+","+centroid+","+offset+","+numsubg+","+not0+","+nameofrows.get(numrows - 1) + ",");
                tablefile_t.write(numrows-1+","+target+","+centroid+","+offset+","+numsubg+","+not0+","+nameofrows.get(numrows - 1) + ",");
                tablefile_c.write(numrows-1+","+target+","+centroid+","+offset+","+numsubg+","+not0+","+nameofrows.get(numrows - 1) + ",");


                //tablefile.write(numrows-1+","+numsubg+","+nameofrows.get(numrows - 1) + ",");

                for (int j = 0; j < numcols-1; j++) {
                    tablefile_u.write(tablefilerow_u[j] + ",");
                    tablefile_s.write(tablefilerow_s[j] + ",");
                    tablefile_t.write(tablefilerow_t[j] + ",");
                    tablefile_c.write(tablefilerow_c[j] + ",");
                }

                tablefile_u.write(tablefilerow_u[numcols-1] + "\n");
                tablefile_s.write(tablefilerow_s[numcols-1] + "\n");
                tablefile_t.write(tablefilerow_t[numcols-1] + "\n");
                tablefile_c.write(tablefilerow_c[numcols-1] + "\n");


                if ((numsubg-not0) != 0) logfile.write("\n***** PROBLEM *****: "+numsubg+" subg but only "+not0+" cols\n");

                //Poi li tolgo
                tablefile_u.flush();
                tablefile_s.flush();
                tablefile_t.flush();
                tablefile_c.flush();

                logfile.flush();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }


        /*
         *
         * THE END
         *
         */

        try {
            logfile.flush();
            logfile.close();
        } catch (IOException e) {
            e.printStackTrace();
        }


        try {
            tablefile_u.flush();
            tablefile_s.flush();
            tablefile_t.flush();
            tablefile_c.flush();

            tablefile_u.close();
            tablefile_s.close();
            tablefile_t.close();
            tablefile_c.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }


    public static void main(String[] args) {
        Driver driver;
        float threshold = 0.4f;

        driver = GraphDatabase.driver("bolt://localhost:7687", AuthTokens.basic("neo4j", "123stella"));
        //LOCAL
        //driver = GraphDatabase.driver("bolt://localhost:7687", AuthTokens.basic("neo4j", "neo4j123Stella"));


        if (args.length == 0) {
            System.out.println("Specify one action! \n Use: java -jar TableCreator <checkstring | table> <clu | shear | tag | combo | all> [threshold 0.x]");
            System.exit(-1);
        }


        if(args[0].compareTo("checkstring") == 0) {
            checkstring(driver);
        } else if(args[0].compareTo("table") == 0) {
            String wlabel = "none";

            if(args[1].compareToIgnoreCase("clu") == 0) wlabel = "weight";
            else if(args[1].compareToIgnoreCase("shear") == 0) wlabel = "shear";
            else if(args[1].compareToIgnoreCase("tag") == 0) wlabel = "tag";
            else if(args[1].compareToIgnoreCase("combo") == 0) wlabel = "combo";
            else if(args[1].compareToIgnoreCase("all") == 0) wlabel = "all";
            else {
                System.out.println("++++ ++++ I cannot determine the kind of weigth you want : is clu, shear, tag, combo or all? ++++ ++++");
                System.exit(-1);
            }

            if(args.length == 3)
                try {
                    threshold =  Float.parseFloat(args[2]);
                } catch (Exception e) {
                    System.out.println("Wrong threshold. I will use 0.4" );
                }

            table(driver,threshold,wlabel);
        } else {
            System.out.println("Wrong action! \n Use: java -jar TableCreator <checkstring | table> <clu | shear | tag | combo | all> [threshold 0.x]");
        }

        driver.close();
        System.out.println("DONE");
        System.exit(0);
    }
}
