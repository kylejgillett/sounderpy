// Splits each documented parameter ("name (type) – description") into a header
// line and a description block, so long parameter lists are easier to scan.
// Pages render normally without it.

document.addEventListener("DOMContentLoaded", () => {
  const DASH = " – ";

  function split(p) {
    const name = p.firstElementChild;
    if (!name || name.tagName !== "STRONG" || p.firstChild !== name) return;

    // Find the dash that separates the type from the description.
    const head = [];
    let node = name.nextSibling, at = -1;
    while (node) {
      if (node.nodeType === Node.TEXT_NODE && (at = node.data.indexOf(DASH)) !== -1) break;
      head.push(node);
      node = node.nextSibling;
    }
    if (!node) return;

    const rest = node.splitText(at);
    rest.data = rest.data.slice(DASH.length);
    if (node.data) head.push(node);

    const type = document.createElement("span");
    type.className = "spy-param-type";
    type.append(...head);
    // Drop the parentheses around the type.
    const first = type.firstChild, last = type.lastChild;
    if (first?.nodeType === Node.TEXT_NODE) first.data = first.data.replace(/^\s*\(/, "");
    if (last?.nodeType === Node.TEXT_NODE) last.data = last.data.replace(/\)\s*$/, "");

    const header = document.createElement("span");
    header.className = "spy-param-head";
    name.classList.add("spy-param-name");
    header.append(name);
    if (type.textContent.trim()) header.append(type);

    const desc = document.createElement("span");
    desc.className = "spy-param-desc";
    for (let n = rest; n; ) {
      const next = n.nextSibling;
      desc.append(n);
      n = next;
    }
    p.replaceChildren(header, desc);
    p.classList.add("spy-param");
  }

  // Folds a "Return type" field into the "Returns" field above it, laid out like a
  // parameter: the type as the header line and the description beneath it.
  function mergeReturns(dl) {
    const label = (dt) => dt.firstChild?.textContent.trim();
    const dts = [...dl.children].filter((el) => el.tagName === "DT");
    const returns = dts.find((dt) => label(dt) === "Returns");
    const rtype = dts.find((dt) => label(dt) === "Return type");
    if (!returns && !rtype) return;

    const head = document.createElement("span");
    head.className = "spy-param-head";
    if (rtype) {
      const type = document.createElement("strong");
      type.className = "spy-param-name spy-return-type";
      type.append(...rtype.nextElementSibling.firstElementChild.childNodes);
      head.append(type);
    }
    const desc = document.createElement("span");
    desc.className = "spy-param-desc";
    if (returns) desc.append(...returns.nextElementSibling.firstElementChild.childNodes);

    const p = document.createElement("p");
    p.className = "spy-param";
    p.append(...[head, desc].filter((el) => el.hasChildNodes()));

    const keep = returns ?? rtype;
    keep.firstChild.textContent = "Returns";
    keep.nextElementSibling.replaceChildren(p);
    if (returns && rtype) {
      rtype.nextElementSibling.remove();
      rtype.remove();
    }
  }

  document.querySelectorAll(".rst-content dl.field-list > dd").forEach((dd) => {
    dd.querySelectorAll(":scope > ul > li > p:first-child, :scope > p:first-child").forEach(split);
  });
  document.querySelectorAll(".rst-content dl.field-list").forEach(mergeReturns);
});
