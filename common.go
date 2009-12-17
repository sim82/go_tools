package common

import (
"os";
"bufio";
)

func IsWhitespace( c byte ) bool {
    return c == ' ' || c == '\t' || c == '\n';
}

func ReadToken( buf string ) (string, string) {
    ptr := 0;
    
    for ; ptr < len(buf) && IsWhitespace(buf[ptr]) ; {
//         fmt.Printf( "skip ws\n" );
        ptr++;
    }
    if ptr == len(buf) {
        return "", buf[ptr:len(buf)];
    }
    start := ptr;
    
    for ; ptr < len(buf) && !IsWhitespace(buf[ptr]) ; {
//         fmt.Printf( "in word\n" );
        ptr++;
    }
    end := ptr;
    
    return buf[start:end], buf[end:len(buf)];
    
}


func LineReader( name string ) chan string {
    out := make( chan string );
    go func() {
        f,_ := os.Open( name, os.O_RDONLY, 0 );
        bio := bufio.NewReader( f );
    
        for ;; {
            s,e := bio.ReadString( '\n' );
        
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
