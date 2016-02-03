
using Docile, Lexicon, Element

# config = Config()
config = Config(md_subheader=:skip, mathjax = true)

index = save("Element.md", Element, config)
save("index.md", Index([index]), config)

# run(`../mkdocs build`)
# mkdocs gh-deploy --clean

