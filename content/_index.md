+++
+++

Correlation matrices are a common way to describe how variables move together. They appear throughout modelling and forecasting: asset returns in finance, genes in biology, sensors in engineering, climate variables in weather models. Each entry is a number between -1 and 1 measuring how strongly two of the variables are related, but the entries are not free to take any values; each correlation must make sense in the context of the others. If A is strongly correlated with B, and B with C, that limits how A and C can correlate. The condition capturing this global consistency is that the matrix be **positive semidefinite**: the correlations must all be achievable by one set of variables at once, not just plausible one pair at a time.

This raises a natural question: given a valid correlation matrix, how far can each entry move before that consistency breaks? Some correlations have plenty of slack while others are pinned down by the rest of the matrix. This site visualizes that slack, to make positive semidefiniteness feel less like an abstract rule and more like a visible geometric constraint. It is built for intuition rather than as a tool for repairing matrices in practice; more principled projection and estimation methods exist for that.

<div class="landing-links">
  <a class="landing-card" href="/explore/">
    <strong>Explore &rarr;</strong>
    <span>Generate or upload a correlation matrix. Click any off-diagonal cell to see its feasible range given the other entries, then drag along the axis to vary its value and watch the other entries' ranges shift in response.</span>
  </a>
  <a class="landing-card" href="/theory/">
    <strong>Theory &rarr;</strong>
    <span>A closed-form derivation of the feasible range, with a partial-correlation interpretation of the bounds.</span>
  </a>
</div>
