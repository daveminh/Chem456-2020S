---
layout: page
title: Tutorials
permalink: /tutorials/
---

<ul id="archive">
{% for tut in site.tutorials %}
      <li class="archiveposturl" style="background: transparent">
        <span><a href="{{ tut.url | prepend: site.baseurl}}">{{ tut.title }}</a></span>
<strong style="font-size:100%; font-family: 'Titillium Web', sans-serif; float:right">
<a title="Download problems (pdf)" href="{{ tut.pdf | prepend: site.baseurl }}"><i class="fas fa-file-pdf"></i></a> 
{% if tut.attachment %}
&nbsp; <a title="Download attachments (zip)" href="{{ tut.attachment | prepend: site.baseurl }}"><i class="fas fa-file-archive"></i></a>
{% endif %}
</strong> 
      </li>
{% endfor %}
</ul>