Submit-block ::= {
  contact {
    contact {
      name name {
        last "",
        first "",
        middle "",
        initials "",
        suffix "",
        title ""
      },
      affil std {
        affil "Laboratory of Marine Biochemical Engineering",
        div "East China University of Science and Technology",
        city "shanghai",
        country "China",
        street "Meilong Road NO.130",
        email "mdna@majorbio.com",
        postal-code "200237"
      }
    }
  },
  cit {
    authors {
      names std {
        {
          name name {
            last "",
            first "",
            middle "",
            initials "",
            suffix "",
            title ""
          }
        }
      },
      affil std {
        affil "Laboratory of Marine Biochemical Engineering",
        div "East China University of Science and Technology",
        city "shanghai",
        country "China",
        street "Meilong Road NO.130",
        postal-code "200237"
      }
    }
  },
  subtype new
}
Seqdesc ::= pub {
  pub {
    gen {
      cit "unpublished",
      authors {
        names std {
          {
            name name {
              last "",
              first "",
              middle "",
              initials "",
              suffix "",
              title ""
            }
          }
        }
      },
      title "whole genome shotgun sequence"
    }
  }
}
Seqdesc ::= user {
  type str "DBLink",
  data {
    {
      label str "BioProject",
      num 1,
      data strs {
        ""
      }
    },
    {
      label str "BioSample",
      num 1,
      data strs {
        ""
      }
    }
  }
}
Seqdesc ::= user {
  type str "Submission",
  data {
    {
      label str "AdditionalComment",
      data str "ALT EMAIL:mdna@majorbio.com"
    }
  }
}
Seqdesc ::= user {
  type str "Submission",
  data {
    {
      label str "AdditionalComment",
      data str "Submission Title:None"
    }
  }
}
