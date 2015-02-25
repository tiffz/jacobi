$(document).ready(function() {
  var matrix = new MatrixWrapper();
  var MAX_OFFB = Math.pow(10, -9);


  // Generate a random symmetric matrix
  $('#generate_matrix').click(function() {
    var number = $('#matrix_size').val() * 1;
    var delimiter = $('#delimiter').val();
    var text = '';
    $('#calc_error').empty();

    var min = $('#max_number').val() * 1;
    var max = $('#min_number').val() * 1;

    matrix = new MatrixWrapper(number);
    matrix.makeRandomSymmetric(min, max);

    var text = matrix.toString(delimiter + ' ');
    $('#matrix_input').val(text);
  });

  $('#calculate_submit').click(function() {
    var error = '';
    var text = $('#matrix_input').val();
    var delimiter = $('#delimiter').val();
    if(text.replace(/\s/g,"") == ""){
      error = "ERROR: You can't have a blank matrix!";
    } else {


      // Parse text into matrix
      matrix = textToMatrixWrapper(text, delimiter);
      $('#matrix_input').val(matrix.toString(delimiter + ' '));
      if (!matrix.isSymmetric()) {
        error = "ERROR: Jacobi's Algorithm only works for symmetrix matrices! ";
      }
    }
    if (!error) {
      advanceSlide('#calculator_tab');
      initializeCalculation(matrix);
    }

    $('#calc_error').html(error);
  });


  $('#simulation_submit').click(function() {
    advanceSlide('#assignment_tab');
    initializeSimulation();
  });

  $('#calculator_tab > .slide:gt(0)').hide();
  $('.calculator_advance').click(function() {
    advanceSlide('#calculator_tab');
  });

  $('#assignment_tab > .slide:gt(0)').hide();
  $('.assignment_advance').click(function() {
    advanceSlide('#assignment_tab');
  });

  function advanceSlide(parent) {
    $(parent + ' > .slide:first').fadeOut()
         .next('.slide').fadeIn()
         .end().appendTo(parent);
  }

  function initializeSimulation() {
    $('#simulation_matrices').html('');
    $('#simulation_navigation').html('');
    var matrices = [];
    var final_matrices = [];
    var m = null;
    var number = 5;
    var min = -10;
    var max = 10;

    var unsortedIterations = [];
    var sortedIterations = [];

    var numMatrices = 10;
    for (var i = 0; i < numMatrices; i++) {
      m = new MatrixWrapper(number);
      m.makeRandomSymmetric(min, max);
      matrices[i] = m;

      // Do sorted iterations
      sortedIterations[i] = 0;
      var offB = m.getOffB();
      var offA = offB;
      var m1 = m;

      var sortedOffB = [offB];
      var unsortedOffB = [offB];

      while(offB > MAX_OFFB) {
        sortedIterations[i]++;
        m1 = m1.iterate();
        offB = m1.getOffB();
        sortedOffB[sortedIterations[i]] = offB;
      }

      // Do unsorted iterations
      unsortedIterations[i] = 0;
      offB = m.getOffB();
      var m2 = m;
      var r = 0;
      var c = r + 1;
      while(offB > MAX_OFFB) {
        unsortedIterations[i]++;
        m2 = m2.iterate(r, c);
        offB = m2.getOffB();
        c++;
        if (c >= m2.cols) {
          r++;
          if (r >= m2.rows - 1) {
            r = 0;
          }
          c = r + 1;
        }

        unsortedOffB[unsortedIterations[i]] = offB;
      }

      var id = (i + 1);
      $('#simulation_navigation').append('<li><a href="#simulation_matrix_' + id + '">Matrix ' + id + '</a></li>');

      $('#simulation_matrices').append('<div class="long block" id="simulation_matrix_' + id + '">'
                                 + '<div class="block"><h2>Original Matrix ' + id
                                 + '</h2>' +  m.prettyPrint() + '</div>'
                                 + '<div class="block"><h2>Final Matrix ' + id
                                 + '</h2>' +  m1.prettyPrint() + '</div>'
                                 + '<div class="long block">Sorted Iterations: ' + sortedIterations[i]
                                 + '; Unsorted Iterations: ' + unsortedIterations[i]
                                 + '</div></div>');
      var sortedId = 'sorted_chart_' + id;
      var unsortedId = 'unsorted_chart_' + id;
      $('#simulation_matrix_' + id).append('<div class="long block"><h2>Plot of OffB over Iterations</h2>'
                                 + '<ul class="simulation_tabs"><li>' 
                                 + '<a href="#' + sortedId + '">With Sorting Step</a></li>'
                                 + '<li><a href="#' + unsortedId + '">Without Sorting Step</a></li></ul>'
                                 + '<canvas id="' + sortedId
                                 + '" class="simulation_chart" width="600px" height="300px"></canvas>'
                                 + '<canvas id="' + unsortedId
                                 + '" class="simulation_chart" width="600px" height="300px"></canvas>'
                                 + '<div class="long small">ln(Off(Bk)) vs k [black] || y = xln(9/10) + ln(Off(A)) [blue]</div>'
                                 + '</div>');
      graph(sortedId, sortedOffB, offA);
      graph(unsortedId, unsortedOffB, offA);


      //alert(sortedIterations[i] + '-' + unsortedIterations[i]);
    }
    $('ul.simulation_tabs').each(function() {
      $(this).simpleTabs();
    });
  }

  function graph(id, offB, offA) {
    var bk = [];
    var k = [];
    var y = [];
    for (var i = 0; i < offB.length; i++) {
      bk[i] = Math.log(offB[i]);
      k[i] = i;
      y[i] = i * Math.log(9 / 10) + Math.log(offA);
    }
    var context = $('#' + id).get(0).getContext("2d");
    var data = {
    labels : k,
    datasets : [
       {
         strokeColor : "#222",
         pointColor : "#222",
         pointStrokeColor : "#eee",
         data : bk
        }, 
       {
         strokeColor : "rgba(151,187,205,1)",
         pointColor : "rgba(151,187,205,1)",
         pointStrokeColor : "#eee",
         data : y
        }
      ]

    }

    var settings = {
                    animation : false, 
                    datasetFill : false, 
                    bezierCurve : false,
                  };

    var chart = new Chart(context).Line(data, settings);
  }

  function initializeCalculation(m) {
    $('#starting_matrix').html(m.prettyPrint());
    $('#iteration_area').empty();
    var iterations = 0;
    var offB = m.getOffB();
    while(offB > MAX_OFFB) {
      iterations++;
      m = m.iterate();
      offB = m.getOffB();
      $('#iteration_area').append('<div class="block"><h2>Iteration ' + iterations
                                 + '</h2>' +  m.prettyPrint() + '<small>OffB: '
                                 + Math.round(offB * 1000000000) / 1000000000 + '</small></div>');
    }
    $("#iteration_area").append('<div id="final_matrix" class="long block">'
             + '<h2>Your Final Matrix</h2>'
             + m.prettyPrint() + '</div>');
    
  }

  $('ul.tabs').each(function() {
    $(this).simpleTabs();
  });
});

/* 
  * jQuery tabs courtesy of Jack Moore
  * www.jacklmoore.com/notes/jquery-tabs/ 
  */
$.fn.simpleTabs = function(){
  // For each set of tabs, we want to keep track of
  // which tab is active and it's associated content
  var $active, $content, $links = $(this).find('a');

  // If the location.hash matches one of the links, use that as the active tab.
  // If no match is found, use the first link as the initial active tab.
  $active = $($links.filter('[href="'+location.hash+'"]')[0] || $links[0]);
  $active.addClass('active');

  $content = $($active[0].hash);

  // Hide the remaining content
  $links.not($active).each(function () {
    $(this.hash).hide();
  });

  // Bind the click event handler
  $(this).on('click', 'a', function(e){
    // Make the old tab inactive.
    $active.removeClass('active');
    $content.hide();

    // Update the variables with the new link and content
    $active = $(this);
    $content = $(this.hash);

    // Make the tab active.
    $active.addClass('active');
    $content.show();

    // Prevent the anchor's default click action
    e.preventDefault();
  });
}