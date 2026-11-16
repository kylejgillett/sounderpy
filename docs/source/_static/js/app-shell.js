// Phones: the sounderpy.app menu button in the docs header (see _templates/layout.html)
document.addEventListener("DOMContentLoaded", () => {
  const root = document.documentElement;
  const toggle = document.querySelector("[data-app-menu-toggle]");
  const menu = document.querySelector("[data-app-menu]");
  if (!toggle || !menu) return;
  const setMenu = (open) => {
    menu.hidden = !open;
    toggle.setAttribute("aria-expanded", String(open));
    toggle.setAttribute("aria-label", open ? "Close menu" : "Open menu");
    root.classList.toggle("app-menu-open", open);
  };
  toggle.addEventListener("click", () => setMenu(menu.hidden));
  menu.addEventListener("click", (e) => { if (e.target.closest("a")) setMenu(false); });
  document.querySelector("[data-app-menu-scrim]")?.addEventListener("click", () => setMenu(false));
  document.addEventListener("keydown", (e) => { if (e.key === "Escape") setMenu(false); });
});
