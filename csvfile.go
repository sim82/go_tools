package main

import (
	"os"
	"fmt"
	"common"
	"strconv"
)

type Column struct {
	//    Type int;
	ColString []string
	ColInt    []int
	ColDouble []float64
}

type CSVFile struct {
	NumLines int
	Cols     []*Column
}

func (cf *CSVFile) GetString(i, j int) string {
	col := cf.Cols[j]
	if col.ColString != nil {
		return col.ColString[i]
	} else {
		panic("not a string column: ", j)
	}
	panic("unreachable")
}

func (cf *CSVFile) GetInt(i, j int) int {
	col := cf.Cols[j]
	if col.ColInt != nil {
		return col.ColInt[i]
	} else {
		panic("not an int column: ", j)
	}
	panic("unreachable")
}


func (cf *CSVFile) GetDouble(i, j int) float64 {
	col := cf.Cols[j]
	if col.ColString != nil {
		return col.ColString[i]
	} else {
		panic("not a double column: ", j)
	}
	panic("unreachable")
}

func ReadCSV(name string, cts string) *CSVFile {
	//     cta := ctaFromString( cts );


	nlines := 0
	for _ = range common.LineReader(name) {
		nlines++
	}

	cf := &CSVFile{
		NumLines: nlines,
		Cols: make([]*Column, len(cts)),
	}

	for i, t := range cts {
		switch t {
		case 'S', 's':
			cf.Cols[i] = &Column{ColString: make([]string, nlines)}
		case 'I', 'i':
			cf.Cols[i] = &Column{ColInt: make([]int, nlines)}
		case 'D', 'd':
			cf.Cols[i] = &Column{ColDouble: make([]float64, nlines)}
		}

	}

	j := 0
	for l := range common.LineReader(name) {
		//lorig := l;
		for i, c := range cf.Cols {
			var token string
			token, l = common.ReadToken(l)

			var e os.Error = nil

			switch {
			case c.ColString != nil:
				cf.Cols[i].ColString[j] = token
			case c.ColInt != nil:
				cf.Cols[i].ColInt[j], e = strconv.Atoi(token)
				if e != nil {
					//fmt.Printf( "line: '%s'\n", string(lorig));
					panic("conversion error from string to int ", i, " ", j, " ", string(token))
				}
			case c.ColDouble != nil:
				cf.Cols[i].ColDouble[j], e = strconv.Atof64(token)
				if e != nil {
					//fmt.Printf( "line: '%s'\n", string(lorig));
					panic("conversion error from string to double ", i, " ", j, " ", string(token))
				}
			}

		}
		j++
	}

	fmt.Printf("nlines: %d\n", nlines)
	return cf
}

func main() {

	cf := ReadCSV("/home/sim/dist_nn/218_nn.txt", "ssiddiiss")

	println(cf.GetString(10, 1))
	_ = cf
}
