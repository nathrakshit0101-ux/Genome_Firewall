CXX      = g++
# C++14 is used for compatibility with MinGW.org GCC 6.3 and modern GCC/Clang.
# -O3 enables aggressive loop unrolling and SIMD vectorization.
CXXFLAGS = -std=c++14 -O3 -pthread \
           -Wall -Wextra \
           -Wno-unused-result -Wno-unused-parameter \
           -Wno-missing-field-initializers \
           -Wno-shift-count-overflow \
           -Iinclude

# Statically link libgcc and libstdc++ so the executable runs standalone on any Windows machine
LDFLAGS  = -pthread -static -static-libgcc -static-libstdc++

SRCDIR   = src
INCDIR   = include
OBJDIR   = obj
BINDIR   = .

SOURCES  = $(SRCDIR)/AMRDatabase.cpp \
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

OBJECTS  = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SOURCES))
TARGET   = genome_firewall

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) $(LDFLAGS) -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJDIR):
	@mkdir -p $(OBJDIR) 2>/dev/null || cmd /c "if not exist $(OBJDIR) mkdir $(OBJDIR)" 2>nul || true

clean:
	@rm -f $(OBJDIR)/*.o $(TARGET).exe $(TARGET) test_suite* 2>/dev/null || cmd /c "if exist $(OBJDIR) del /Q $(OBJDIR)\*.o & if exist $(TARGET).exe del /Q $(TARGET).exe" 2>nul || true

test: $(TARGET)
	@echo [*] Running Verification Test on CARD Top Markers and Sample Genome...
	./$(TARGET) --input data/sample_real.fasta --db data/markers_card_top.csv

unit-test: tests/test_genome_firewall.cpp $(filter-out $(OBJDIR)/main.o $(OBJDIR)/REPLShell.o, $(OBJECTS))
	$(CXX) $(CXXFLAGS) tests/test_genome_firewall.cpp $(filter-out $(OBJDIR)/main.o $(OBJDIR)/REPLShell.o, $(OBJECTS)) $(LDFLAGS) -o test_suite
	./test_suite
	@rm -f test_suite.exe test_suite 2>/dev/null || true

.PHONY: all clean test unit-test
