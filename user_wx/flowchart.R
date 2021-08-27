library(DiagrammeR)
grViz("digraph flowchart {
      # node definitions with substituted label text
      node [fontname = Helvetica, shape = rectangle]        
      tab1 [label = '@@1']
      tab2 [label = '@@2']
      tab3 [label = '@@3']
      tab4 [label = '@@4']
      tab5 [label = '@@5']
      # edge definitions with the node IDs
      tab5 -> tab1
      tab1 -> tab3 
      tab2 -> tab3 
      tab3 -> tab4;
      }

      [1]: 'dataset fixed'
      [2]: 'treatment and control input variable'
      [3]: 'model from our range of models created'
      [4]: 'results'
      [5]: 'raw counts'
      ")


grViz("digraph flowchart {
      # node definitions with substituted label text
      node [fontname = Helvetica, shape = rectangle]        
      tab1 [label = '@@1']
      tab2 [label = '@@2']
      tab3 [label = '@@3']
      tab4 [label = '@@4']
      tab5 [label = '@@5']
      # edge definitions with the node IDs
      tab5 -> tab1
      tab2 -> tab1 
      tab2 -> tab3 
      tab1 -> tab3
      tab3 -> tab4;
      }

      [1]: 'dataset variable'
      [2]: 'treatment and control input variable'
      [3]: 'model from our range of models created'
      [4]: 'results'
      [5]: 'raw counts'
      ")

grViz("digraph flowchart {
      # node definitions with substituted label text
      node [fontname = Helvetica, shape = rectangle]        
      tab1 [label = '@@1']
      tab2 [label = '@@2']
      tab3 [label = '@@3']
      tab4 [label = '@@4']
      tab5 [label = '@@5']
      tab6 [label = '@@6']
      # edge definitions with the node IDs

      tab2 -> tab1 
      tab1 -> tab4
      tab4 -> tab3 
      tab2 -> tab3
      tab5 -> tab1
      tab3 -> tab6;
      }

      [1]: 'dataset variable'
      [2]: 'treatment and control input variable'
      [3]: 'PCA'
      [4]: 'gene signature variable'
      [5]: 'raw counts'
      [6]: 'same signature PCA components diff for different treatment and controls'
      ")
grViz("digraph flowchart {
      # node definitions with substituted label text
      node [fontname = Helvetica, shape = rectangle]        
      tab1 [label = '@@1']
      tab2 [label = '@@2']
      tab3 [label = '@@3']

      tab5 [label = '@@5']

      tab7 [label = '@@7']
      # edge definitions with the node IDs

      tab2 -> tab1 
      tab5 -> tab1
      tab1 -> tab7
      tab3 -> tab7
      tab2 -> tab7;
      }

      [1]: 'dataset variable'
      [2]: 'treatment and control input variable'
      [3]: 'mediator'

      [5]: 'raw counts'

      [7]: 'mediation model'
      ")
