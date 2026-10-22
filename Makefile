CC       = gcc
CXX      = g++
CFLAGS   = -O3 -Wall -Wextra -Iinclude
CXXFLAGS = -std=c++14 -O3 -pthread \
           -Wall -Wextra \
           -Wno-unused-result -Wno-unused-parameter \
           -Wno-missing-field-initializers \
           -Wno-shift-count-overflow \
           -Iinclude

LDFLAGS  = -pthread -static -static-libgcc -static-libstdc++

SRCDIR   = src
INCDIR   = include
OBJDIR   = obj
BINDIR   = .

C_SOURCES   = $(SRCDIR)/nucleotide_bitpack.c
CPP_SOURCES = $(SRCDIR)/AMRDatabase.cpp \
              $(SRCDIR)/FastaMmapParser.cpp \
              $(SRCDIR)/AMRBloomFilter.cpp \
              $(SRCDIR)/TrieArena.cpp \
              $(SRCDIR)/AMRTrie.cpp \
              $(SRCDIR)/LogisticRegressor.cpp \
              $(SRCDIR)/PredictorEngine.cpp \
              $(SRCDIR)/ScoringEngine.cpp \
              $(SRCDIR)/ProteomeTranslator.cpp \
              $(SRCDIR)/REPLShell.cpp \
              $(SRCDIR)/main.cpp

OBJECTS  = $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(C_SOURCES)) \
           $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(CPP_SOURCES))

TARGET   = genome_firewall

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) $(LDFLAGS) -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.c | $(OBJDIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJDIR):
	@mkdir -p $(OBJDIR) 2>/dev/null || cmd /c "if not exist $(OBJDIR) mkdir $(OBJDIR)" 2>nul || true

clean:
	@rm -f $(OBJDIR)/*.o $(TARGET).exe $(TARGET) test_suite* fasta_benchmark* card_extractor* 2>/dev/null || cmd /c "if exist $(OBJDIR) del /Q $(OBJDIR)\*.o & if exist $(TARGET).exe del /Q $(TARGET).exe" 2>nul || true

test: $(TARGET)
	@echo [*] Running Verification Test on CARD Top Markers and Sample Genome...
	./$(TARGET) --input data/sample_real.fasta --db data/markers_card_top.csv

unit-test: tests/test_genome_firewall.cpp $(filter-out $(OBJDIR)/main.o $(OBJDIR)/REPLShell.o, $(OBJECTS))
	$(CXX) $(CXXFLAGS) tests/test_genome_firewall.cpp $(filter-out $(OBJDIR)/main.o $(OBJDIR)/REPLShell.o, $(OBJECTS)) $(LDFLAGS) -o test_suite
	./test_suite
	@rm -f test_suite.exe test_suite 2>/dev/null || true

benchmark: tools/fasta_benchmark.cpp $(filter-out $(OBJDIR)/main.o $(OBJDIR)/REPLShell.o, $(OBJECTS))
	$(CXX) $(CXXFLAGS) tools/fasta_benchmark.cpp $(filter-out $(OBJDIR)/main.o $(OBJDIR)/REPLShell.o, $(OBJECTS)) $(LDFLAGS) -o fasta_benchmark
	./fasta_benchmark
	@rm -f fasta_benchmark.exe fasta_benchmark 2>/dev/null || true

extract-tool: tools/card_dataset_extractor.cpp
	$(CXX) $(CXXFLAGS) tools/card_dataset_extractor.cpp -o card_extractor
	./card_extractor
	@rm -f card_extractor.exe card_extractor 2>/dev/null || true

.PHONY: all clean test unit-test benchmark extract-tool
