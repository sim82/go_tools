package main

import ( "fmt"; "bufio"; "os"; "strconv";  )


type MultipleAlignment struct {
    names []string;
    data [][]byte;
    
    index map[string]int;
}

type MultipleAlignmentIf interface {
    buildIndex();
}

func LineReader( name string ) chan []byte {
    
    f,_ := os.Open( name, os.O_RDONLY, 0 );
    bio := bufio.NewReader( f );
    
    out := make( chan []byte );
    
    go func() {
        for ;; {
            s,e := bio.ReadSlice( '\n' );
        
            if e == nil {
                out <- s;
            } else {
                close(out);
                break;
            }
        }
        
        f.Close();
        return;
    }();
    
    return out;
    
}

func isWhitespace( c byte ) bool {
    return c == ' ' || c == '\t' || c == '\n';
}

func readToken( buf []byte ) ([]byte, []byte) {
    ptr := 0;
    
    for ; ptr < len(buf) && isWhitespace(buf[ptr]) ; {
//         fmt.Printf( "skip ws\n" );
        ptr++;
    }
    if ptr == len(buf) {
        return []byte{}, buf[ptr:len(buf)];
    }
    start := ptr;
    
    for ; ptr < len(buf) && !isWhitespace(buf[ptr]) ; {
//         fmt.Printf( "in word\n" );
        ptr++;
    }
    end := ptr;
    
    return buf[start:end], buf[end:len(buf)];
    
}

func (self *MultipleAlignment) buildIndex() {
    self.index = map[string]int {};
    
    for i := 0; i < len(self.names); i++ {
        self.index[self.names[i]] = i;
    }
}

func (self *MultipleAlignment) dataByName( name string ) []byte {
    if self.index == nil {
        panic( "no index" );
    }
    
    i := self.index[name];
    return self.data[i];
}

func main() {
    
    ma := ReadPhylip( "/space/dist_subseq_alignments/1604_0000_200_60" );
    ma.buildIndex();
    
    d := ma.dataByName( "CbnHyd17_01" );
    
    fmt.Printf( "%s\nlen: %d\n", string(d), len(d) );
}

func ReadPhylip( name string ) *MultipleAlignment {
    first := true;
    
    var nTaxa, seqLen int;
    
    var names []string;
    var datii [][]byte;
    
    var i int;
    for line := range LineReader(name) {
        if first {
            var e os.Error;
            
            var t []byte;
            t,line = readToken( line );
           
            nTaxa,e = strconv.Atoi( string(t) );
            
            if e != nil {
                panic( fmt.Sprintf( "cannot interpret nTaxa field: %s\n", string(t) ));
            }
            //            fmt.Printf( "rest: %s\n", line );
            t,line = readToken( line );
            
            //            fmt.Printf( "token: '%s'\n", t );
            seqLen,e = strconv.Atoi( string(t) );
            if e != nil {
                panic( fmt.Sprintf( "cannot interpret seqLen field: %s\n", string(t) ));
            }
            
            names = make( []string, nTaxa );
            datii = make( [][]byte, nTaxa );
            
            first = false;
        } else {
            var tName []byte;
            tName,line = readToken( line );
            
          
            var start int;
            for start = 0; start < len(line) && isWhitespace(line[start]) ; {
                start++;
            }
            
            var end int;
            for end = len(line) - 1; end > 0 && isWhitespace(line[end]) ; {
                end--;
            }
            
            if start >= end {
                panic( "cannot read sequence data\n" );
            }
            
            
            data := line[start:end+1];
            
            if len(data) != seqLen {
                panic( "len(data) != seqLen" );
            }
            
            names[i] = string(tName);
            datii[i] = data;
            i++;
            fmt.Printf( "taxon: %s %d\n", tName, len(data));
        }
    }
    
    fmt.Printf( "nTaxa: %d len: %d\n", nTaxa, seqLen );
 
    return &MultipleAlignment { names, datii, nil };
}
