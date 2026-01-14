#include <tinyxml2.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <regex>
#include <cstdlib>
#include <algorithm>

using namespace tinyxml2;
using namespace std;

/* -------------------------------
   Unimod mapping
-------------------------------- */
std::map<std::string, std::string> UNIMOD_MAP = {

   /* =========================
      Phosphorylation
      ========================= */
   {"Phosphoserine",        "UNIMOD:21"},
   {"Phosphothreonine",     "UNIMOD:21"},
   {"Phosphotyrosine",      "UNIMOD:21"},

   /* =========================
      Acetylation / Acylation
      ========================= */
   {"Acetylation",         "UNIMOD:1"},
   {"N-acetylation",       "UNIMOD:1"},
   {"N-acetylserine",      "UNIMOD:1"},
   {"N-acetylalanine",     "UNIMOD:1"},
   {"N-acetylmethionine",  "UNIMOD:1"},
   {"N-acetylglycine",     "UNIMOD:1"},
   {"N-acetylproline",     "UNIMOD:1"},
   {"N-acetylthreonine",   "UNIMOD:1"},
   {"N6-acetyllysine",     "UNIMOD:1"},
   {"Propionylation",      "UNIMOD:58"},
   {"Butyrylation",        "UNIMOD:128"},
   {"Crotonylation",       "UNIMOD:136"},
   {"Succinylation",       "UNIMOD:64"},
   {"Malonylation",        "UNIMOD:747"},
   {"Glutarylation",       "UNIMOD:121"},

   /* =========================
      Methylation
      ========================= */
   {"Methylation",          "UNIMOD:34"},
   {"Dimethylation",        "UNIMOD:36"},
   {"Trimethylation",       "UNIMOD:37"},

   /* =========================
      Hydroxylation
      ========================= */
   {"Hydroxylation",        "UNIMOD:35"},
   {"4-hydroxyproline",     "UNIMOD:35"},
   {"5-hydroxylysine",      "UNIMOD:35"},

   /* =========================
      Glycosylation (generic)
      ========================= */
   {"N-linked glycosylation", "UNIMOD:43"},
   {"O-linked glycosylation", "UNIMOD:43"},
   {"Glycosylation",          "UNIMOD:43"},
   {"HexNAc",                 "UNIMOD:43"},

   /* =========================
      Lipidation
      ========================= */
   {"Myristoylation",       "UNIMOD:211"},
   {"Palmitoylation",       "UNIMOD:214"},
   {"Prenylation",          "UNIMOD:235"},
   {"Farnesylation",        "UNIMOD:134"},
   {"Geranylgeranylation",  "UNIMOD:206"},

   /* =========================
      Sulfation / Nitrosylation
      ========================= */
   {"Sulfation",            "UNIMOD:40"},
   {"Nitrosylation",        "UNIMOD:36"},

   /* =========================
      Ubiquitin-like modifiers
      ========================= */
   {"Ubiquitination",       "UNIMOD:121"},
   {"Sumoylation",          "UNIMOD:314"},
   {"NEDDylation",          "UNIMOD:384"},

   /* =========================
      ADP-ribosylation
      ========================= */
   {"ADP-ribosylation",     "UNIMOD:384"},

   /* =========================
      Terminal modifications
      ========================= */
   {"Amidation",            "UNIMOD:2"},
   {"Formylation",          "UNIMOD:122"},

   /* =========================
      Redox / cysteine biology
      ========================= */
   {"Disulfide bond",       "UNIMOD:7"},
   {"Glutathionylation",    "UNIMOD:55"}
};

/* ============================================================ */
unordered_map<string, size_t> PTM_COUNTS;
unordered_map<string, size_t> VARIANT_SKIPPED;
unordered_map<string, size_t> VARIANT_COMPLEX;

vector<string> OUTPUT_SIMPLE;


/* ============================================================
   Helpers
   ============================================================ */
static bool is_valid_aa(char c)
{
   return (c >= 'A' && c <= 'Z' && c != 'B' && c != 'J' && c != 'O' && c != 'U' && c != 'X' && c != 'Z');
}

// Safe string getter - returns empty string if pointer is null
static string safe_string(const char* str)
{
   return str ? string(str) : string();
}

/* ============================================================
   Parse modified residues
   ============================================================ */
vector<pair<int, string>> parse_modres(XMLElement* entry, bool strict)
{
   vector<pair<int, string>> mods;

   for (auto* feat = entry->FirstChildElement("feature");
      feat; feat = feat->NextSiblingElement("feature"))
   {
      // Safe attribute access with default value
      string feat_type = safe_string(feat->Attribute("type"));
      if (feat_type != "modified residue")
         continue;

      const char* desc = feat->Attribute("description");
      auto* locElem = feat->FirstChildElement("location");
      auto* posElem = locElem ? locElem->FirstChildElement("position") : nullptr;

      // Skip if description or position is missing
      if (!desc || !posElem) continue;

      string d(desc);
      PTM_COUNTS[d]++;

      auto it = UNIMOD_MAP.find(d);
      if (it == UNIMOD_MAP.end())
      {
         if (strict)
         {
            cerr << "ERROR: Unmapped PTM: " << d << endl;
            exit(1);
         }
         continue;
      }

      mods.emplace_back(posElem->IntAttribute("position"), it->second);
   }
   return mods;
}

/* ============================================================
   Parse VariantSimple (Comet-safe)
   ============================================================ */
vector<string> parse_variants(XMLElement* entry)
{
   vector<string> variants;
   regex simple_sub(R"(([A-Z])\s*->\s*([A-Z]))");

   for (auto* feat = entry->FirstChildElement("feature");
      feat; feat = feat->NextSiblingElement("feature"))
   {
      // Safe attribute access with default value
      string feat_type = safe_string(feat->Attribute("type"));
      if (feat_type != "sequence variant")
         continue;

      const char* desc = feat->Attribute("description");
      auto* locElem = feat->FirstChildElement("location");
      auto* posElem = locElem ? locElem->FirstChildElement("position") : nullptr;

      // Get original and variation elements
      auto* origElem = feat->FirstChildElement("original");
      auto* varElem = feat->FirstChildElement("variation");

      // Skip if position is missing
      if (!posElem)
      {
         VARIANT_SKIPPED["complex_location"]++;
         continue;
      }

      // Check if SGRP variant (skip these)
      if (desc && string(desc).find("SGRP") != string::npos)
      {
         VARIANT_SKIPPED["SGRP"]++;
         continue;
      }

      char ref = '\0';
      char alt = '\0';

      // Try to parse from <original> and <variation> elements first
      if (origElem && varElem)
      {
         const char* origText = origElem->GetText();
         const char* varText = varElem->GetText();

         if (origText && varText && strlen(origText) == 1 && strlen(varText) == 1)
         {
            ref = origText[0];
            alt = varText[0];
         }
      }

      // If that didn't work, try parsing from description
      if (ref == '\0' && desc)
      {
         string d(desc);
         smatch m;
         if (regex_search(d, m, simple_sub))
         {
            ref = m[1].str()[0];
            alt = m[2].str()[0];
         }
      }

      // If we still don't have valid ref/alt, skip
      if (ref == '\0' || alt == '\0')
      {
         VARIANT_SKIPPED["non_simple"]++;
         continue;
      }

      if (!is_valid_aa(ref) || !is_valid_aa(alt))
      {
         VARIANT_SKIPPED["invalid_aa"]++;
         continue;
      }

      int pos = posElem->IntAttribute("position");
      variants.push_back("(" + to_string(pos) + "|" + ref + "|" + alt + ")");
   }

   return variants;
}

/* ============================================================
   Parse VariantComplex
   Handles: insertions, deletions, multi-residue substitutions,
            and splice variants
   Format: Each variant is (StartPos|EndPos|Sequence)
   - Deletions: (pos1|pos2|) - empty sequence
   - Insertions: (pos|pos|INSERTED_SEQ)
   - Substitutions: (pos1|pos2|VARIANT_SEQ)
   ============================================================ */
vector<string> parse_complex_variants(XMLElement* entry)
{
   vector<string> simple_vars;
   vector<string> complex_vars;

   for (auto* feat = entry->FirstChildElement("feature"); feat; feat = feat->NextSiblingElement("feature"))
   {
      string feat_type = safe_string(feat->Attribute("type"));

      if (feat_type != "sequence variant" && feat_type != "splice variant" && feat_type != "mutagenesis site")
      {
         continue;
      }

      const char* desc = feat->Attribute("description");

      auto* locElem = feat->FirstChildElement("location");

      if (!locElem)
      {
         VARIANT_SKIPPED["complex_location"]++;
         continue;
      }

      auto* origElem = feat->FirstChildElement("original");
      auto* varElem = feat->FirstChildElement("variation");

      const char* origText = origElem ? origElem->GetText() : nullptr;
      const char* varText = varElem ? varElem->GetText() : nullptr;

      if (!origText || !varText)
      {
         VARIANT_SKIPPED["non_simple"]++;
         continue;
      }

      auto* posElem = locElem->FirstChildElement("position");
      auto* beginElem = locElem->FirstChildElement("begin");
      auto* endElem = locElem->FirstChildElement("end");

      int start_pos = 0;
      int end_pos = 0;

      if (posElem)
      {
         start_pos = end_pos = posElem->IntAttribute("position");
      }
      else if (beginElem && endElem)
      {
         start_pos = beginElem->IntAttribute("position");
         end_pos = endElem->IntAttribute("position");
      }
      else if (beginElem || endElem)
      {
         start_pos = beginElem ? beginElem->IntAttribute("position") : endElem->IntAttribute("position");
         end_pos = endElem ? endElem->IntAttribute("position") : start_pos;
      }
      else
      {
         VARIANT_SKIPPED["complex_location"]++;
         continue;
      }

      // Determine if this should be VariantSimple
      if (safe_string(origText).length() == 1 && safe_string(varText).length() == 1 &&
         is_valid_aa(origText[0]) && is_valid_aa(varText[0]) &&
         start_pos == end_pos)
      {
         string simple_entry = "(" + to_string(start_pos) + "|" + safe_string(origText) + "|" + safe_string(varText) + ")";
         simple_vars.push_back(simple_entry);
         continue; // Skip processing as VariantComplex
      }

      // If not simple, classify as VariantComplex
      string complex_entry = "(" + to_string(start_pos) + "|" + to_string(end_pos) + "|" + safe_string(varText) + ")";
      complex_vars.push_back(complex_entry);

      if (feat_type == "mutagenesis site")
      {
         VARIANT_COMPLEX["mutagenesis"]++;
      }
      else if (feat_type == "sequence variant")
      {
         VARIANT_COMPLEX["sequence_variant"]++;
      }
   }

   // Pass simple_vars to OUTPUT_SIMPLE
   OUTPUT_SIMPLE = simple_vars;

   return complex_vars;
}

/* ============================================================
   Main
   ============================================================ */
int main(int argc, char* argv[])
{
   if (argc < 3)
   {
      cerr << "Usage: " << argv[0] << " input.xml output.peff [options]\n";
      cerr << "Options:\n";
      cerr << "  --strict              Exit on unmapped PTMs (default: skip)\n";
      cerr << "  --no-ptms             Disable PTM processing (default: enabled)\n";
      cerr << "  --variant-simple      Enable VariantSimple processing (default: disabled)\n";
      cerr << "  --variant-complex     Enable VariantComplex processing (default: disabled)\n";
      return 1;
   }

   // Parse command-line options
   bool strict = false;
   bool enable_ptms = true;
   bool enable_variant_simple = false;
   bool enable_variant_complex = false;

   for (int i = 3; i < argc; i++)
   {
      string arg(argv[i]);
      if (arg == "--strict")
      {
         strict = true;
      }
      else if (arg == "--no-ptms")
      {
         enable_ptms = false;
      }
      else if (arg == "--variant-simple")
      {
         enable_variant_simple = true;
      }
      else if (arg == "--variant-complex")
      {
         enable_variant_complex = true;
      }
      else
      {
         cerr << "Unknown option: " << arg << endl;
         return 1;
      }
   }

   XMLDocument doc;
   if (doc.LoadFile(argv[1]) != XML_SUCCESS)
   {
      cerr << "Failed to load XML\n";
      return 1;
   }

   ofstream peff(argv[2]);
   if (!peff.is_open())
   {
      cerr << "Failed to open output file: " << argv[2] << endl;
      return 1;
   }

   peff << "# PEFF 1.0\n";
   peff << "# Database=UniProt\n";
   peff << "# VariantSimple=" << (enable_variant_simple || enable_variant_complex ? "true" : "false") << "\n";
   peff << "# VariantComplex=" << (enable_variant_complex ? "true" : "false") << "\n";
   peff << "# ModResUnimod=" << (enable_ptms ? "true" : "false") << "\n";

   XMLElement* root = doc.RootElement();
   if (!root)
   {
      cerr << "Invalid XML: no root element\n";
      return 1;
   }

   size_t entries_processed = 0;
   size_t entries_skipped = 0;

   for (auto* entry = root->FirstChildElement("entry"); entry; entry = entry->NextSiblingElement("entry"))
   {
      auto* accElem = entry->FirstChildElement("accession");
      auto* seqElem = entry->FirstChildElement("sequence");

      if (!accElem || !seqElem)
      {
         cerr << "Entry missing ";
         if (!accElem)
         {
            cerr << "<accession> ";
         }
         if (!seqElem)
         {
            cerr << "<sequence> ";
         }
         cerr << "- processing features within this entry regardless" << endl;
      }

      string acc = accElem ? accElem->GetText() : "UNKNOWN_ACCESSION";
      string seq = seqElem ? seqElem->GetText() : "UNKNOWN_SEQUENCE";

      seq.erase(remove(seq.begin(), seq.end(), '\n'), seq.end());
      seq.erase(remove(seq.begin(), seq.end(), ' '), seq.end());

      string db_type = "tr"; // Default to TrEMBL
      const char* dataset = entry->Attribute("dataset");
      if (dataset && string(dataset) == "Swiss-Prot")
      {
         db_type = "sp";
      }

      string entry_name = safe_string(entry->FirstChildElement("name") ? entry->FirstChildElement("name")->GetText() : "");
      string organism = safe_string(entry->FirstChildElement("organism") ? entry->FirstChildElement("organism")->GetText() : "");

      // Clear OUTPUT_SIMPLE for this entry
      OUTPUT_SIMPLE.clear();

      // Parse features based on enabled options
      vector<pair<int, string>> mods;
      vector<string> vars;
      vector<string> complex_vars;

      if (enable_ptms)
      {
         mods = parse_modres(entry, strict);
      }

      if (enable_variant_simple)
      {
         vars = parse_variants(entry);
      }

      if (enable_variant_complex)
      {
         complex_vars = parse_complex_variants(entry);
      }

      peff << ">" << db_type << "|" << acc << "|" << entry_name;
      if (!organism.empty())
      {
         peff << " OS=" << organism;
      }

      // Write VariantSimple
      // When VariantComplex is enabled, prefer OUTPUT_SIMPLE from parse_complex_variants (for compatibility)
      // Otherwise, use vars from parse_variants if VariantSimple is enabled
      vector<string> simple_output;
      if (enable_variant_complex && !OUTPUT_SIMPLE.empty())
      {
         simple_output = OUTPUT_SIMPLE;
      }
      else if (enable_variant_simple)
      {
         simple_output = vars;
      }
      
      if (!simple_output.empty())
      {
         peff << " \\VariantSimple=";
         for (size_t i = 0; i < simple_output.size(); ++i)
         {
            peff << simple_output[i];
         }
      }

      // Write VariantComplex
      if (!complex_vars.empty())
      {
         peff << " \\VariantComplex=";
         for (size_t i = 0; i < complex_vars.size(); ++i)
         {
            peff << complex_vars[i];
         }
      }

      if (!mods.empty())
      {
         peff << " \\ModResUnimod=";
         for (size_t i = 0; i < mods.size(); i++)
         {
            peff << "(" << mods[i].first << "|" << mods[i].second << ")";
         }
      }

      peff << "\n";
      if (!seq.empty() && seq != "UNKNOWN_SEQUENCE")
      {
         peff << seq << "\n";
      }

      entries_processed++;
   }

   peff.close();

   cerr << "Done. Processed " << entries_processed << " entries, skipped " << entries_skipped << " entries.\n";

   return 0;
}
