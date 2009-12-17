/*
 * Copyright (C) 2009 Simon A. Berger
 * 
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 2 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 */


/**
 * A go implementation of a reader for phylogenetic trees in the newick format
 * (http://evolution.genetics.washington.edu/phylip/newicktree.html)
 *
 * For now it should be possible to read *unrooted* standard newick files (= the
 * top level node must have three children!). It can also read branch-labels and
 * node labels which will be interpreted as branch support values if possible
 * 
 * There is also code for a basic tree writer.
 * See func main for usage information.
 */

package main

import (
"fmt";
"strconv";
"io";
"strings";
"os";
)

type AData struct {
    label string;
    isTip bool;
    support float;
}

type LN struct {
    next *LN;
    back *LN;
    
    backLen float;
    backLabel string;
    backSupport float;
    data *AData;
}


func CreateLN() (*LN){
    data := &AData{};
    var ns [3]LN;
    ns[0].next = &ns[1];
    ns[1].next = &ns[2];
    ns[2].next = &ns[0];
    ns[0].data = data;
    ns[1].data = data;
    ns[2].data = data;
    
    return &ns[0];
    
    //     data := &AData{};
//     
//     n := &LN{&LN{&LN{nil,nil,0,"",0,data}, nil, 0, "", 0, data}, nil, 0, "", 0, data};
//     n.next.next.next = n;
//     
//     return n;


}

func CreateLeaf( name string ) (*LN) {
    n := CreateLN();
    
    n.data.label = name;
    n.data.isTip = true;
    
    return n;
}

type ParserInput interface {
  // Get the character at an arbitrary position in the
  // input source.
  CharAt(i int) uint8;

  // Get the number of characters in the input source.
  Size() int;
  
  Substring( start, end int ) (string)
}

func IsDigit( c uint8 ) (bool) {
    return c >= '0' && c <= '9';
}

func IsSpace( c uint8 ) (bool) {
    return c == 32 || c == 9;
}

func skipWhitespace( pi ParserInput, pos int) ( int ) {
    for ;;pos++ {
        if !IsSpace(pi.CharAt( pos )) {
            break;
        }
    }
   
    return pos;
}

func findEndOfBranch( pi ParserInput, pos int ) (int) {
    termchar := ":,)";
    
    for ;; pos++ {
        for c := range termchar {
            if termchar[c] == pi.CharAt(pos) {
                return pos;
            }
        }
    }
    
    panic( "unreachable" );
}


func isFloatChar(c uint8) (bool){
    return (c>='0' && c <= '9') || c == '.' || c == 'e' || c == 'E' || c == '-';
}

func findFloat( pi ParserInput, pos int ) (int){
    
    
    for ; isFloatChar(pi.CharAt(pos));  {
        pos++;
    }
    //fmt.Printf( "float len: %d\n", len );
    return pos;
}
        
func parseBranchLength( pi ParserInput, pos int ) (float, int) {
    pos = skipWhitespace( pi, pos );

    // expect + consume ':'
    if pi.CharAt(pos) != ':' {
        //throw new RuntimeException("parse error: parseBranchLength expects ':' at " + ptr);
        return 0.0, 0;
    } else {
        
        pos++;
        
        pos = skipWhitespace(pi, pos);
        
        flen := findFloat( pi, pos);
        if flen == pos {
            fmt.Printf("error: missing float number at %d\n", pos );
        }
        
        blen,_ := strconv.Atof(pi.Substring(pos, flen));
        pos = flen;
        return blen,pos;
       
    }
    
    panic("unreachable");
}

func findNext( pi ParserInput, pos int, c uint8  ) (int) {
    for ;;pos++ {
        if pi.CharAt( pos ) == c {
            break;
        }
        
    }
    return pos;
}


func parseBranchLabel( pi ParserInput, pos int ) (string, int) {
    var label string;
    if( pi.CharAt(pos) == '[' ) {
        lstart := pos;
        pos++;
        
        
        lend := findNext( pi, pos, ']' );
        
        pos = lend+1;
        
        // labels have the form [blablabla], so the label content starts at lstart + 1
        
        if lend - (lstart+1) <= 0  {
            //printLocation();
            fmt.Printf("bad branch label: %s\n", pi.Substring(lstart, lend) );
        }
        
        label = pi.Substring( lstart + 1, lend);
        
    } else {
        label = "";
    }
    
    return label, pos;
}
    
    
func twiddle( n1 *LN, n2 *LN, len float, label string, support float ) {
    n1.back = n2;
    n2.back = n1;
    n1.backLen = len;
    n2.backLen = len;
    n1.backLabel = label;
    n2.backLabel = label;
    n1.backSupport = support;
    n2.backSupport = support;
}
    
func parseInnerNode( pi ParserInput, pos int ) (*LN, int) {

    //fmt.Printf( "parse inner node\n" );
    
    pos = skipWhitespace( pi, pos );
    
    if pi.CharAt( pos ) != '(' {
        fmt.Printf( "error: expected '('\n" );
    }
    // consume opening '('
    pos++;
    
    var nl *LN;
    nl, pos = parseNode( pi, pos );
   
    var bll float;
    bll, pos = parseBranchLength( pi, pos );
 
//     fmt.Printf( "parsing next node %d\n", pos );
    
    var labell string;
    labell, pos = parseBranchLabel( pi, pos );
    
    pos = skipWhitespace( pi, pos );
    
    
    // expect + consume ','
    if (pi.CharAt( pos ) != ',') {
        
        fmt.Printf("parse error: parseInnerNode expects ',' at %d\n", pos );
    }
    pos++;

    // parse right node + branch length
    var nr *LN;
    nr, pos = parseNode( pi, pos );
    
    var blr float;
    blr, pos = parseBranchLength( pi, pos );
    
    var labelr string;
    labelr, pos = parseBranchLabel( pi, pos );
    
    pos = skipWhitespace( pi, pos );
    
    if( pi.CharAt( pos ) == ',') {
        // second comma found: three child nodes == pseudo root
        pos++;
        
        var nx *LN;
        nx, pos = parseNode( pi, pos );
        
        var blx float;
        blx, pos= parseBranchLength( pi, pos );

        var labelx string;
        labelx, pos = parseBranchLabel( pi, pos );
       
        pos = skipWhitespace( pi, pos );
        
        if( pi.CharAt( pos ) != ')' ) {
            //            printLocation();
            fmt.Printf("parse error: parseInnerNode (at root) expects ') at %d\n", pos );
        }
        pos++;
        pos = skipWhitespace( pi, pos );
        
        n := CreateLN();
        
        twiddle( nl, n.next, bll, labell, n.data.support );
        twiddle( nr, n.next.next, blr, labelr, n.data.support );
        twiddle( nx, n, blx, labelx, n.data.support );
        
  //      fmt.Printf( "twiddle root\n" );
        
        //                      System.out.printf( "root: %f %f %f\n", nl.data.getSupport(), nr.data.getSupport(), nx.data.getSupport() );
        //                      System.exit(0);
        return n, pos;
    } else if pi.CharAt( pos ) == ')' {
        // the stuff between the closing '(' and the ':' of the branch length
        // is interpreted as node-label. If the node label corresponds to a float value
        // it is interpreted as branch support (or node support as a rooted-trees-only-please biologist would say)
        
        pos++;
        lend := findEndOfBranch(pi, pos);
        
        nodeLabel := pi.Substring(pos, lend);
        pos = lend;
                
    //    fmt.Printf( "pos: %d %s\n", pos, nodeLabel );
        isDigit := true;
        for i := range nodeLabel {
            isDigit = isDigit && IsDigit(nodeLabel[i]);
            
            if( i == 0 ) {
                isDigit = isDigit && (nodeLabel[i] != '0');
            }
        }
        
        
        var support float;
        if( isDigit ) {
            support,_ = strconv.Atof(nodeLabel);
        } else {
            support = -1;
        }
        
        n := CreateLN();
        n.data.support = support;
        //n.data.setNodeLabel(nodeLabel);
        
        twiddle( nl, n.next, bll, labell, support );
        twiddle( nr, n.next.next, blr, labelr, support );
        
        return n,pos;
      
        
    } else {
        fmt.Printf("parse error: parseInnerNode expects ')'or ',' at %d got: %c\n", pos, pi.CharAt( pos ));
        return nil,pos;
    }
    
    panic( "unreachable" );
}

func parseLeaf( pi ParserInput, pos int ) (*LN, int) {
    pos = skipWhitespace( pi, pos );
    
    brLen := findEndOfBranch( pi, pos );
    
    ld := pi.Substring(pos,brLen);
    pos = brLen;
    
    node := CreateLeaf( ld );
    
    //fmt.Printf( "leaf parsed: %s\n", ld );
    return node, pos;
}

func parseNode( pi ParserInput, pos int ) (*LN, int) {
    pos = skipWhitespace( pi, pos );
  
    var n *LN;
    
    if pi.CharAt( pos ) == '(' {
        n, pos = parseInnerNode( pi, pos );
    } else {
        n, pos = parseLeaf( pi, pos );
    }
    return n,pos;
}

func Parse( pi ParserInput, pos int ) ( *LN, int ) {
    pos = skipWhitespace( pi, pos );
    
    var n *LN;
    n, pos = parseNode( pi, pos );
    
    return n, pos;
    
}
type StringPI string;

func (self StringPI) CharAt( i int ) (uint8) {
//    fmt.Printf( "CharAt: %d: %d %c\n", i, self[i], self[i]);
    return (self)[i];
}

func (self StringPI) Size() (int) {
    return len(self);
}

func (self StringPI) Substring( start, end  int ) (string) {
    return string(self)[start:end];
}


func printTreeInner( node *LN, w io.Writer, root bool ) {
    if( node.data.isTip ) {
        //s.printf("%s:%G", node.data.getTipName(), node.backLen);
        w.Write( strings.Bytes(fmt.Sprintf("%s:%8.20f", node.data.label, node.backLen)));
    } else {
        w.Write( []byte{'('} );
        printTreeInner( node.next.back, w, false);
        w.Write([]byte{','});
        printTreeInner(node.next.next.back, w, false);

        if( root ) {
            w.Write([]byte{','});
            printTreeInner( node.back, w, false);
            w.Write( []byte{')',';'} );
        } else {
            //s.printf("):%G", node.backLen );
            w.Write( strings.Bytes(fmt.Sprintf("):%8.20f", node.backLen )));
        }
    }
}

func PrintTree( node *LN, w io.Writer ) {
    
    
    //node = LN.getTowardsTree(node);
    if( node.data.isTip ) {
        if( node.back != nil ) {
            node = node.back;
        } else if( node.next.back != nil ) {
            node = node.next.back;
        } else if( node.next.next.back != nil ) {
            node = node.next.next.back;
        } else {
            panic( "can not print single unlinked node");
        }
        
        if( node.data.isTip ) {
            panic( "could not find non-tip node for writing the three (this is a braindead limitation of this tree printer!)");
        }
    }
    printTreeInner( node, w, true );
}
    


func main() {
    data,_ := io.ReadFile("/space/raxml/VINCENT/RAxML_bipartitions.1604.BEST.WITH" );
    
    n,_ := Parse( StringPI( string( data )), 0);
    //n,_ := Parse( StringPI(" (bla:.1,bla2:2.1,(bla32:1.2,bla34:0.7)100:3.4);"), 0 );
    
    PrintTree( n, os.Stdout );
    fmt.Println();
    PrintTree( n.back, os.Stdout );
    fmt.Println();
}
