EXEC = tree_parser csvfile
LIB = common.6 tree.6

tree_parser.6: tree.6

%.6: %.go
	6g -I. $(<)

%: %.6 $(LIB)
	6l -L. -o $(@) $(<)

all: $(EXEC)

clean:
	rm -rf $(LIB)
	rm -rf $(EXEC)